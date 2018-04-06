#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_incompressible.h"

#include <Python.h>  // has to be the first included header
#include <memory>
#include <petscksp.h>
#include <petscsys.h>
#include <cmath>

#include "semt/Semt.h"
#include "semt/Shortcuts.h"
#include "control/types.h"
#include "spatial_discretization/finite_element_method/solid_mechanics/elasticity_tensor.h"

namespace SpatialDiscretization
{
  
// general implementation for solid mechanics
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  Mesh::isDeformable<typename BasisOnMeshType::Mesh>,
  Equation::isIncompressible<Term>
>::
setStiffnessMatrix()
{
  // get pointer to mesh object
  std::shared_ptr<BasisOnMeshType> mesh = this->data_.mesh();
  
  const int D = BasisOnMeshType::dim();  // = 3
  const int nDofs = mesh->nDofs();
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  const int nElements = mesh->nElements();
  const int nUnknowsPerElement = nDofsPerElement*D;    // 3 directions for displacements per dof
  
  // define shortcuts for quadrature
  typedef Quadrature::TensorProduct<D,QuadratureType> QuadratureDD;
  
  // define type to hold evaluations of integrand for stiffness matrix
  typedef MathUtility::Matrix<nUnknowsPerElement,nUnknowsPerElement> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureType::numberEvaluations()
          > EvaluationsArrayType;     // evaluations[nGP^D](nUnknows,nUnknows)
  
  // setup arrays used for integration
  std::array<std::array<double,D>, QuadratureType::numberEvaluations()> samplingPoints = QuadratureType::samplingPoints();
  EvaluationsArrayType evaluationsArray;
  
  typedef std::array<Vec3,BasisOnMeshType::dim()> Tensor2;
  PetscErrorCode ierr;
  Mat &tangentStiffnessMatrix = this->data_.tangentStiffnessMatrix();
  
  // initialize values to zero
  int cntr = 1;
  // loop over elements 
  for (element_no_t elementNo = 0; elementNo < mesh->nElements(); elementNo++)
  {
    auto dofNo = mesh->getElementDofNos(elementNo);
    
    for (int aDof = 0; aDof < nDofsPerElement; aDof++)
    {
      for (int aComponent = 0; aComponent < D; aComponent++)
      {
        for (int bDof = 0; bDof < nDofsPerElement; bDof++)
        {
          for (int bComponent = 0; bComponent < D; bComponent++)
          {
            dof_no_t matrixRowIndex = dofNo[aDof]*D + aComponent;
            dof_no_t matrixColumnIndex = dofNo[bDof]*D + bComponent;
        
            VLOG(3) << " initialize tangentStiffnessMatrix entry ( " << matrixRowIndex << "," << matrixColumnIndex << ") (no. " << cntr++ << ")";
            ierr = MatSetValue(tangentStiffnessMatrix, matrixRowIndex, matrixColumnIndex, 0, INSERT_VALUES); CHKERRV(ierr);
          }
        }
      }
    }
  }
  
  // loop over elements 
  for (int elementNo = 0; elementNo < nElements; elementNo++)
  { 
    // get geometry field of reference configuration
    std::array<Vec3,nDofsPerElement> geometryReferenceValues;
    this->data_.geometryReference().getElementValues(elementNo, geometryReferenceValues);
    
    // get displacement field
    std::array<Vec3,nDofsPerElement> displacementValues;
    this->data_.displacements().getElementValues(elementNo, displacementValues);
    
    // loop over integration points (e.g. gauss points) for displacement field
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // get parameter values of current sampling point
      std::array<double,D> xi = samplingPoints[samplingPointIndex];
      
      // compute the 3xD jacobian of the parameter space to world space mapping
      Tensor2 jacobianMaterial = BasisOnMeshType::computeJacobian(geometryReferenceValues, xi);
      double jacobianDeterminant;
      Tensor2 inverseJacobianMaterial = MathUtility::computeInverse(jacobianMaterial, jacobianDeterminant);
      // jacobianMaterial[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx
      // inverseJacobianMaterial[columnIdx][rowIdx] = dxi_rowIdx/dX_columnIdx because of inverse function theorem
      
      // F
      Tensor2 deformationGradient = this->computeDeformationGradient(displacementValues, inverseJacobianMaterial, xi);
      double deformationGradientDeterminant = MathUtility::computeDeterminant(deformationGradient);  // J
      
      Tensor2 rightCauchyGreen = this->computeRightCauchyGreenTensor(deformationGradient);  // C = F^T*F
      
      double rightCauchyGreenDeterminant;   // J^2
      Tensor2 inverseRightCauchyGreen = MathUtility::computeSymmetricInverse(rightCauchyGreen, rightCauchyGreenDeterminant);  // C^-1
      
      // invariants
      std::array<double,3> invariants = this->computeInvariants(rightCauchyGreen, rightCauchyGreenDeterminant);  // I_1, I_2, I_3
      std::array<double,2> reducedInvariants = this->computeReducedInvariants(invariants, deformationGradientDeterminant); // Ibar_1, Ibar_2
  
      // artifical pressure p
      double artificialPressureTilde;
      const double artificialPressure = this->computeArtificialPressure(deformationGradientDeterminant, artificialPressureTilde);
      
      // Pk2 stress tensor S = S_vol + S_iso (p.234)
      std::array<Vec3,3> fictitiousPK2Stress;
      std::array<Vec3,3> pk2StressIsochoric;
      Tensor2 PK2Stress = this->computePK2Stress(artificialPressure, rightCauchyGreen, inverseRightCauchyGreen, reducedInvariants, deformationGradientDeterminant, fictitiousPK2Stress, pk2StressIsochoric);      
      // elasticity tensor C_{ijkl}
      ElasticityTensor elasticity = this->computeElasticityTensor(artificialPressure, artificialPressureTilde, rightCauchyGreen, inverseRightCauchyGreen, fictitiousPK2Stress, pk2StressIsochoric, deformationGradientDeterminant);
      
      std::array<Vec3,nDofsPerElement> gradPhi = mesh->getGradPhi(xi);
      // (column-major storage) gradPhi[M][a] = dphi_M / dxi_a
      // gradPhi[column][row] = gradPhi[dofIndex][i] = dphi_dofIndex/dxi_i, columnIdx = dofIndex, rowIdx = which direction
      
      // loop over pairs of basis functions and evaluate integrand at xi
      for (int aDof = 0; aDof < nDofsPerElement; aDof++)           // index over dofs, each dof has D components, L in derivation
      {
        for (int aComponent = 0; aComponent < D; aComponent++)     // lower-case a in derivation, index over displacement components
        {
          // compute index of degree of freedom and component (matrix row index)
          const int i = nDofsPerElement*aDof + aComponent;
         
          for (int bDof = 0; bDof < nDofsPerElement; bDof++)       // index M in derivation
          {
            for (int bComponent = 0; bComponent < D; bComponent++)     // lower-case b in derivation, index over displacement components
            {
         
              // compute index of degree of freedom and component (index of entry in vector of unknows / matrix column index)
              const int j = nDofsPerElement*bDof + bComponent;
             
              // compute integrand[i][j]
              double integrand = 0.0;
              for (int dInternal = 0; dInternal < D; dInternal++)     // capital D in derivation
              {
                for (int bInternal = 0; bInternal < D; bInternal++)     // capital B in derivation
                {
                  // ----------------------------
                  // compute derivatives of phi
                  const double dphiL_dXB = 0.0;
                  const double dphiM_dXD = 0.0;
                  
                  // helper index k for multiplication with inverse Jacobian
                  for (int k = 0; k < D; k++)
                  {
                    // compute dphiL/dXB from dphiL/dxik and dxik/dXB
                    const double dphiL_dxik = gradPhi[aDof][k];    // dphi_L/dxik
                    const double dxik_dXB = inverseJacobianMaterial[bInternal][k];  // inverseJacobianMaterial[B][k] = J_kB = dxi_k/dX_B
                    
                    dphiL_dXB += dphiL_dxik * dxik_dXB;
                    
                    // compute dphiM/dXD from dphiM/dxik and dxik/dXD
                    const double dphiM_dxik = gradPhi[bDof][k];    // dphi_M/dxik
                    const double dxik_dXD = inverseJacobianMaterial[dInternal][k];  // inverseJacobianMaterial[D][k] = J_kD = dxi_k/dX_D
                    
                    dphiM_dXD += dphiM_dxik * dxik_dXD;
                  }
                 
                 
                  // ----------------------------
                  // compute ktilde_abBD
                  const int delta_ab = (aComponent == bComponent? 1 : 0);
                  const double S_BD = PK2Stress[dInternal][bInternal];
                  
                  // compute ffC = FaA * FbC * C_ABCD
                  const double ffC = 0.0;
                  for (int aInternal = 0; aInternal < D; aInternal++)     // capital A in derivation
                  {
                    for (int cInternal = 0; cInternal < D; cInternal++)    // capital C in derivation
                    {
                      const double FaA = deformationGradient[aInternal][aComponent];
                      const double FbC = deformationGradient[cInternal][bComponent];
                      const double cc = elasticity.getEntry(aInternal, bInternal, cInternal, dInternal);
                      ffC += FaA * FbC * cc;
                    }
                  }
                  
                  const double ktilde_abBD = delta_ab*S_BD + ffC;
                 
                  // k_ij = int dphi_La/dX_B * dphi_Mb/dX_D * ktilde_abBD
                  integrand += dphiL_dXB * ktilde_abBD * dphiM_dXD;
                }
              }
              evaluationsArray[samplingPointIndex](i,j) = integrand * fabs(jacobianDeterminant);
            }  // b
          }  // M 
        }  // a
      }  // L
    }  // function evaluations
    
    // integrate all values for the (i,j) dof pairs at once
    EvaluationsType integratedValues = QuadratureDD::computeIntegral(evaluationsArray);
    
    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNo = mesh->getElementDofNos(elementNo);
    
    // set entries in stiffness matrix
    // loop over indices of unknows ((aDof,aComponent), (bDof,bComponent)), i.e. (i,j)
    for (int aDof = 0; aDof < nDofsPerElement; aDof++)
    {
      for (int aComponent = 0; aComponent < D; aComponent++)
      {
        // compute index of degree of freedom and component (matrix row index)
        const int i = nDofsPerElement*aDof + aComponent;
       
        for (int bDof = 0; bDof < nDofsPerElement; bDof++)
        {
          for (int bComponent = 0; bComponent < D; bComponent++)
          {
            // compute index of degree of freedom and component (index of entry in vector of unknows / matrix column index)
            const int j = nDofsPerElement*bDof + bComponent;
                
            // integrate value and set entry in stiffness matrix
            double integratedValue = integratedValues(i,j);
            
            // Compute indices in stiffness matrix. For each dof there are D=3 values for the 3 displacement components. 
            // Therefore the nDofsPerElement number is not the number of unknows.
            dof_no_t matrixRowIndex = dofNo[aDof]*D + aComponent;
            dof_no_t matrixColumnIndex = dofNo[bDof]*D + bComponent;
            
            VLOG(2) << "  unknown pair (("<<aDof<<","<<aComponent<<"),("<<bDof<<","<<bComponent<<")) = (" <<i<<","<<j<<"), dofs (" << dofNo[aDof] << ","<<dofNo[bDof]<<")";
            VLOG(2) << "      matrix indices ("<<matrixRowIndex<<","<<matrixColumnIndex<<"), integrated value: "<<integratedValue;
            
            ierr = MatSetValue(tangentStiffnessMatrix, matrixRowIndex, matrixColumnIndex, integratedValue, ADD_VALUES); CHKERRV(ierr);
            
          }  // bComponent
        }  // bDof       
      }  // aComponent
    }  // aDof
  }  // elementNo
  
  // because this is used in nonlinear solver context, assembly has to be called here, not via data->finalAssembly
  MatAssemblyBegin(tangentStiffnessMatrix,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(tangentStiffnessMatrix,MAT_FINAL_ASSEMBLY);
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  Mesh::isDeformable<typename BasisOnMeshType::Mesh>,
  Equation::isIncompressible<Term>
>::
setDisplacements(Vec &x)
{
  // copy values of x to displacements
  VecCopy(x,this->data_.displacements());
  
  // set entries of Dirichlet BCs to specified values
  VecSetValues(this->data_.displacements(), dirichletIndices_.data(), dirichletValues_.data(), INSERT_VALUES);
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
Vec &FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  Mesh::isDeformable<typename BasisOnMeshType::Mesh>,
  Equation::isIncompressible<Term>
>::
tangentStiffnessMatrix()
{
  return this->data_.tangentStiffnessMatrix();
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  Mesh::isDeformable<typename BasisOnMeshType::Mesh>,
  Equation::isIncompressible<Term>
>::
computeInternalVirtualWork(Vec &resultVec)
{
  // get pointer to mesh object
  std::shared_ptr<BasisOnMeshType> mesh = this->data_.mesh();
  
  const int D = BasisOnMeshType::dim();  // = 3
  const int nDofs = mesh->nDofs();
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  const int nElements = mesh->nElements();
  const int nUnknowsPerElement = nDofsPerElement*D;    // 3 directions for displacements per dof
  
  // define shortcuts for quadrature
  typedef Quadrature::TensorProduct<D,QuadratureType> QuadratureDD;
  
  // define type to hold evaluations of integrand for stiffness matrix
  typedef std::array<double, nUnknowsPerElement> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureType::numberEvaluations()
          > EvaluationsArrayType;     // evaluations[nGP^D](nUnknows,nUnknows)
  
  // setup arrays used for integration
  std::array<std::array<double,D>, QuadratureType::numberEvaluations()> samplingPoints = QuadratureType::samplingPoints();
  EvaluationsArrayType evaluationsArray;
  
  typedef std::array<Vec3,BasisOnMeshType::dim()> Tensor2;
  PetscErrorCode ierr;
  
  // get memory from PETSc where to store result (not used)
  //const PetscScalar *result;
  //VecGetArray(resultVec, &result);
  
  // loop over elements 
  for (int elementNo = 0; elementNo < nElements; elementNo++)
  { 
    // get geometry field of reference configuration
    std::array<Vec3,nDofsPerElement> geometryReferenceValues;
    this->data_.geometryReference().getElementValues(elementNo, geometryReferenceValues);
    
    // get displacement field
    std::array<Vec3,nDofsPerElement> displacementValues;
    this->data_.displacements().getElementValues(elementNo, displacementValues);
    
    // loop over integration points (e.g. gauss points) for displacement field
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // get parameter values of current sampling point
      std::array<double,D> xi = samplingPoints[samplingPointIndex];
      
      // compute the 3xD jacobian of the parameter space to world space mapping
      Tensor2 jacobianMaterial = BasisOnMeshType::computeJacobian(geometryReferenceValues, xi);
      double jacobianDeterminant;
      Tensor2 inverseJacobianMaterial = MathUtility::computeInverse(jacobianMaterial, jacobianDeterminant);
      // jacobianMaterial[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx
      // inverseJacobianMaterial[columnIdx][rowIdx] = dxi_rowIdx/dX_columnIdx because of inverse function theorem
      
      // F
      Tensor2 deformationGradient = this->computeDeformationGradient(displacementValues, inverseJacobianMaterial, xi);
      double deformationGradientDeterminant = MathUtility::computeDeterminant(deformationGradient);  // J
      
      Tensor2 rightCauchyGreen = this->computeRightCauchyGreenTensor(deformationGradient);  // C = F^T*F
      
      double rightCauchyGreenDeterminant;   // J^2
      Tensor2 inverseRightCauchyGreen = MathUtility::computeSymmetricInverse(rightCauchyGreen, rightCauchyGreenDeterminant);  // C^-1
      
      // invariants
      std::array<double,3> invariants = this->computeInvariants(rightCauchyGreen, rightCauchyGreenDeterminant);  // I_1, I_2, I_3
      std::array<double,2> reducedInvariants = this->computeReducedInvariants(invariants, deformationGradientDeterminant); // Ibar_1, Ibar_2
  
      // artifical pressure p
      double artificialPressureTilde;
      const double artificialPressure = this->computeArtificialPressure(deformationGradientDeterminant, artificialPressureTilde);
      
      // Pk2 stress tensor S = S_vol + S_iso (p.234)
      std::array<Vec3,3> fictitiousPK2Stress;
      std::array<Vec3,3> pk2StressIsochoric;
      Tensor2 PK2Stress = this->computePK2Stress(artificialPressure, rightCauchyGreen, inverseRightCauchyGreen, reducedInvariants, deformationGradientDeterminant, fictitiousPK2Stress, pk2StressIsochoric);      
      
      std::array<Vec3,nDofsPerElement> gradPhi = mesh->getGradPhi(xi);
      // (column-major storage) gradPhi[M][a] = dphi_M / dxi_a
      // gradPhi[column][row] = gradPhi[dofIndex][i] = dphi_dofIndex/dxi_i, columnIdx = dofIndex, rowIdx = which direction
      
      // loop basis functions and evaluate integrand at xi
      for (int aDof = 0; aDof < nDofsPerElement; aDof++)           // index over dofs, each dof has D components, L in derivation
      {
        for (int aComponent = 0; aComponent < D; aComponent++)     // lower-case a in derivation, index over displacement components
        {
          // compute index of degree of freedom and component (result vector index)
          const int i = nDofsPerElement*aDof + aComponent;
         
          // compute result[i]
          double integrand = 0.0;
          for (int aInternal = 0; aInternal < D; aInternal++)     // capital A in derivation
          {
            for (int bInternal = 0; bInternal < D; bInternal++)     // capital B in derivation
            {
              const double faB = deformationGradient[bInternal][aComponent];
              const double faA = deformationGradient[aInternal][aComponent];
              
              // ----------------------------
              // compute derivatives of phi
              const double dphiL_dXA = 0.0;
              const double dphiL_dXB = 0.0;
              
              // helper index k for multiplication with inverse Jacobian
              for (int k = 0; k < D; k++)
              {
                // compute dphiL/dXA from dphiL/dxik and dxik/dXA
                const double dphiL_dxik = gradPhi[aDof][k];    // dphi_L/dxik
                const double dxik_dXA = inverseJacobianMaterial[aInternal][k];  // inverseJacobianMaterial[A][k] = J_kA = dxi_k/dX_A
                
                dphiL_dXA += dphiL_dxik * dxik_dXA;
                
                // compute dphiL/dXB from dphiL/dxik and dxik/dXB
                const double dxik_dXB = inverseJacobianMaterial[bInternal][k];  // inverseJacobianMaterial[B][k] = J_kB = dxi_k/dX_B
                
                dphiL_dXB += dphiL_dxik * dxik_dXB;
              }
              
              integrand += PK2Stress[bInternal][aInternal] * (faB*dphiL_dXA + faA*dphiL_dXB);
                         
            }  // B, bInternal
          }  // A, aInternal
        }  // a
      }  // L
    }  // function evaluations
    
    // integrate all values for result vector entries at once
    EvaluationsType integratedValues = QuadratureDD::computeIntegral(evaluationsArray);
    
    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNo = mesh->getElementDofNos(elementNo);
    
    // add entries in result vector
    // loop over indices of unknows (aDof,aComponent)
    for (int aDof = 0; aDof < nDofsPerElement; aDof++)
    {
      for (int aComponent = 0; aComponent < D; aComponent++)
      {
        // compute index of degree of freedom and component (matrix row index)
        const int i = nDofsPerElement*aDof + aComponent;
   
        // integrate value and set entry in stiffness matrix
        double integratedValue = integratedValues[i];
        
        // Compute indices in stiffness matrix. For each dof there are D=3 values for the 3 displacement components. 
        // Therefore the nDofsPerElement number is not the number of unknows.
        dof_no_t resultVectorIndex = dofNo[aDof]*D + aComponent;
        
        VLOG(2) << "  result vector "<<aDof<<","<<aComponent<<" = " <<i<<", dof " << dofNo[aDof];
        VLOG(2) << "      vector index "<<resultVectorIndex<<", integrated value: "<<integratedValue;
            
        ierr = VecSetValue(resultVec, resultVectorIndex, integratedValue, ADD_VALUES); CHKERRV(ierr);
        
      }  // aComponent
    }  // aDof
  }  // elementNo
  
  VecAssemblyBegin(resultVec);
  VecAssemblyEnd(resultVec);
  // return memory acces of result vector back to PETSc (not used)
  //VecRestoreArray(resultVec, &result);
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  Mesh::isDeformable<typename BasisOnMeshType::Mesh>,
  Equation::isIncompressible<Term>
>::
initialize()
{
  this->data_.initialize();
  this->setStiffnessMatrix();
  this->setRightHandSide();
  this->data_.finalAssembly();
  initializeBoundaryConditions();
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  Mesh::isDeformable<typename BasisOnMeshType::Mesh>,
  Equation::isIncompressible<Term>
>::
initializeBoundaryConditions()
{
  LOG(TRACE)<<"initializeBoundaryConditions";
 
  dof_no_t nUnknowns = this->data_.nUnknowns();
  Vec &rightHandSide = this->data_.externalVirtualEnergy();
  
  PetscErrorCode ierr;
  
  // traverse Dirichlet boundary conditions
  
  // get the first dirichlet boundary condition from the list
  std::pair<node_no_t, double> boundaryCondition 
    = PythonUtility::getOptionDictBegin<dof_no_t, double>(this->specificSettings_, "DirichletBoundaryCondition");
  
  // loop over Dirichlet boundary conditions
  for (; !PythonUtility::getOptionDictEnd(this->specificSettings_, "DirichletBoundaryCondition"); 
       PythonUtility::getOptionDictNext<dof_no_t, double>(this->specificSettings_, "DirichletBoundaryCondition", boundaryCondition))
  {
    dof_no_t boundaryConditionUnknownsIndex = boundaryCondition.first;
    double boundaryConditionValue = boundaryCondition.second;
    
    if (boundaryConditionUnknownsIndex < 0)
      continue;
    
    if (boundaryConditionUnknownsIndex > nUnknowns) 
    {
      LOG(WARNING) << "Boundary condition specified for degree of freedom no. "<<boundaryConditionUnknownsIndex
       <<", but scenario has only "<<nUnknowns<<" degrees of freedom.";
       continue;
    }
    
    // get the rhs entry of the prescribed entry
    double rhsValue;
    ierr = VecGetValues(rightHandSide, 1, &boundaryConditionUnknownsIndex, &rhsValue); CHKERRV(ierr);

    // store values
    dirichletIndices_.push_back(boundaryConditionUnknownsIndex);
    dirichletValues_.push_back(boundaryConditionValue);
    rhsValues_.push_back(rhsValue);
  }
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  Mesh::isDeformable<typename BasisOnMeshType::Mesh>,
  Equation::isIncompressible<Term>
>::
applyDirichletBoundaryConditionsInNonlinearFunction(Vec &f)
{
  // overwrite the values in f that correspond to entries for which Dirichlet BC are set with rhs values, such that nonlinear equation is satisfied for them
  VecSetValues(f, dirichletIndices_.size(), dirichletIndices_.data(), rhsValues_.data(), INSERT_VALUES);
}
  
};    // namespace