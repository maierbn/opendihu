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
#include "basis_on_mesh/00_basis_on_mesh_base_dim.h"
#include "utility/python_utility.h"

namespace SpatialDiscretization
{
  
// general implementation for solid mechanics penalty 
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
setStiffnessMatrix(Mat stiffnessMatrix)
{
  // get pointer to mesh object
  std::shared_ptr<BasisOnMeshType> mesh = this->data_.mesh();
  
  const int D = BasisOnMeshType::dim();  // = 3
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  const int nElements = mesh->nElements();
  const int nUnknowsPerElement = nDofsPerElement*D;    // 3 directions for displacements per dof
  
  // define shortcuts for quadrature
  typedef Quadrature::TensorProduct<D,QuadratureType> QuadratureDD;
  
  // define type to hold evaluations of integrand for stiffness matrix
  typedef MathUtility::Matrix<nUnknowsPerElement,nUnknowsPerElement> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureDD::numberEvaluations()
          > EvaluationsArrayType;     // evaluations[nGP^D](nUnknows,nUnknows)
  
  // setup arrays used for integration
  std::array<VecD<D>, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();
  EvaluationsArrayType evaluationsArray;
  
  PetscErrorCode ierr;
  Mat &tangentStiffnessMatrix = (stiffnessMatrix == PETSC_NULL ? this->data_.tangentStiffnessMatrix() : stiffnessMatrix);
  
  LOG(TRACE) << "setStiffnessMatrix - compute tangent stiffness matrix";
  LOG(DEBUG) << "nUnknowsPerElement: " << nUnknowsPerElement<<", n evaluations for quadrature: " << QuadratureDD::numberEvaluations();
  
  // initialize values to zero
  if (!tangentStiffnessMatrixInitialized_)
  {
    // there are no previous non-zero entries, so loop over all later needed entries and manually set to zero
   
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
          
              VLOG(3) << " initialize tangentStiffnessMatrix (("<<aDof<<","<<aComponent<<"),("<<bDof<<","<<bComponent<<")), dofs (" << dofNo[aDof] << ","<<dofNo[bDof]<<"), entry ( " << matrixRowIndex << "," << matrixColumnIndex << ") (no. " << cntr++ << ")";
              ierr = MatSetValue(tangentStiffnessMatrix, matrixRowIndex, matrixColumnIndex, 0.0, INSERT_VALUES); CHKERRV(ierr);
            }
          }
        }
      }
    }
  }
  else
  {
    // zero all previous non-zero entries, keeps matrix structure
    MatZeroEntries(tangentStiffnessMatrix);
    
    VLOG(3) << " zero tangentStiffnessMatrix";              
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
      VecD<D> xi = samplingPoints[samplingPointIndex];
      
      // compute the 3xD jacobian of the parameter space to world space mapping
      Tensor2<D> jacobianMaterial = MathUtility::transformToDxD<D,D>(BasisOnMeshType::computeJacobian(geometryReferenceValues, xi));
      double jacobianDeterminant;
      Tensor2<D> inverseJacobianMaterial = MathUtility::computeInverse<D>(jacobianMaterial, jacobianDeterminant);
      // jacobianMaterial[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx
      // inverseJacobianMaterial[columnIdx][rowIdx] = dxi_rowIdx/dX_columnIdx because of inverse function theorem
      
      // F
      Tensor2<D> deformationGradient = this->computeDeformationGradient(displacementValues, inverseJacobianMaterial, xi);
      double deformationGradientDeterminant = MathUtility::computeDeterminant<D>(deformationGradient);  // J
      
      Tensor2<D> rightCauchyGreen = this->computeRightCauchyGreenTensor(deformationGradient);  // C = F^T*F
      
      double rightCauchyGreenDeterminant;   // J^2
      Tensor2<D> inverseRightCauchyGreen = MathUtility::computeSymmetricInverse<D>(rightCauchyGreen, rightCauchyGreenDeterminant);  // C^-1
      
      // invariants
      std::array<double,3> invariants = this->computeInvariants(rightCauchyGreen, rightCauchyGreenDeterminant);  // I_1, I_2, I_3
      std::array<double,2> reducedInvariants = this->computeReducedInvariants(invariants, deformationGradientDeterminant); // Ibar_1, Ibar_2
  
      // artifical pressure p
      double artificialPressureTilde;
      const double artificialPressure = this->computeArtificialPressure(deformationGradientDeterminant, artificialPressureTilde);
      
      // Pk2 stress tensor S = S_vol + S_iso (p.234)
      Tensor2<D> fictitiousPK2Stress;
      Tensor2<D> pk2StressIsochoric;
      Tensor2<D> PK2Stress = this->computePK2Stress(artificialPressure, rightCauchyGreen, inverseRightCauchyGreen, reducedInvariants, deformationGradientDeterminant, fictitiousPK2Stress, pk2StressIsochoric);      
      // elasticity tensor C_{ijkl}
      ElasticityTensor elasticity = this->computeElasticityTensor(artificialPressure, artificialPressureTilde, rightCauchyGreen, inverseRightCauchyGreen, fictitiousPK2Stress, pk2StressIsochoric, deformationGradientDeterminant, reducedInvariants);

      std::array<VecD<D>,nDofsPerElement> gradPhi = mesh->getGradPhi(xi);
      // (column-major storage) gradPhi[M][a] = dphi_M / dxi_a
      // gradPhi[column][row] = gradPhi[dofIndex][i] = dphi_dofIndex/dxi_i, columnIdx = dofIndex, rowIdx = which direction
      
      VLOG(2) << "geometryReferenceValues: " << geometryReferenceValues;
      VLOG(2) << "Jacobian: J_phi=" << jacobianMaterial;
      VLOG(2) << "jacobianDeterminant: J=" << jacobianDeterminant;
      VLOG(2) << "inverseJacobianMaterial: J_phi^-1=" << inverseJacobianMaterial;
      VLOG(2) << "deformationGradient: F=" << deformationGradient;
      VLOG(2) << "deformationGradientDeterminant: det F=" << deformationGradientDeterminant;
      VLOG(2) << "rightCauchyGreen: C=" << rightCauchyGreen;
      VLOG(2) << "rightCauchyGreenDeterminant: det C=" << rightCauchyGreenDeterminant;
      VLOG(2) << "inverseRightCauchyGreen: C^-1=" << inverseRightCauchyGreen;
      VLOG(2) << "invariants: I1,I2,I3: " << invariants;
      VLOG(2) << "reducedInvariants: Ibar1, Ibar2: " << reducedInvariants;
      VLOG(2) << "artificialPressure: p=" << artificialPressure << ", artificialPressureTilde: pTilde=" << artificialPressureTilde;
      VLOG(2) << "fictitiousPK2Stress: Sbar=" << fictitiousPK2Stress;
      VLOG(2) << "pk2StressIsochoric: S_iso=" << pk2StressIsochoric;
      VLOG(2) << "PK2Stress: S=" << PK2Stress;
      VLOG(2) << "gradPhi: " << gradPhi;
      VLOG(2) << "elasticity: C=" << elasticity;
      
      // loop over pairs of basis functions and evaluate integrand at xi
      for (int aDof = 0; aDof < nDofsPerElement; aDof++)           // index over dofs, each dof has D components, L in derivation
      {
        for (int aComponent = 0; aComponent < D; aComponent++)     // lower-case a in derivation, index over displacement components
        {
          // compute index of degree of freedom and component (matrix row index)
          const int i = D*aDof + aComponent;
         
          for (int bDof = 0; bDof < nDofsPerElement; bDof++)       // index M in derivation
          {
            for (int bComponent = 0; bComponent < D; bComponent++)     // lower-case b in derivation, index over displacement components
            {
         
              // compute index of degree of freedom and component (index of entry in vector of unknows / matrix column index)
              const int j = D*bDof + bComponent;
             
              // compute integrand[i][j]
              double integrand = 0.0;
              for (int dInternal = 0; dInternal < D; dInternal++)     // capital D in derivation
              {
                for (int bInternal = 0; bInternal < D; bInternal++)     // capital B in derivation
                {
                  // ----------------------------
                  // compute derivatives of phi
                  double dphiL_dXB = 0.0;
                  double dphiM_dXD = 0.0;
                  
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
                  double ffC = 0.0;
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
              
              //LOG(DEBUG) << "(L,a),(M,b)=(" << aDof << ","<<aComponent<<"),("<<bDof<<","<<bComponent<<"), (i,j)=("<<i<<","<<j<<"), integrand=" << integrand;
              
              // store integrand in evaluations array
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
        const int i = D*aDof + aComponent;
       
        for (int bDof = 0; bDof < nDofsPerElement; bDof++)
        {
          for (int bComponent = 0; bComponent < D; bComponent++)
          {
            // compute index of degree of freedom and component (index of entry in vector of unknows / matrix column index)
            const int j = D*bDof + bComponent;
                
            // integrate value and set entry in stiffness matrix
            double integratedValue = integratedValues(i,j);
            
            // Compute indices in stiffness matrix. For each dof there are D=3 values for the 3 displacement components. 
            // Therefore the nDofsPerElement number is not the number of unknows.
            dof_no_t matrixRowIndex = dofNo[aDof]*D + aComponent;
            dof_no_t matrixColumnIndex = dofNo[bDof]*D + bComponent;
            
            //VLOG(2) << "  unknown pair (("<<aDof<<","<<aComponent<<"),("<<bDof<<","<<bComponent<<")) = (" <<i<<","<<j<<"), dofs (" << dofNo[aDof] << ","<<dofNo[bDof]<<")";
            //VLOG(2) << "      matrix indices ("<<matrixRowIndex<<","<<matrixColumnIndex<<"), integrated value: "<<integratedValue;
            
            ierr = MatSetValue(tangentStiffnessMatrix, matrixRowIndex, matrixColumnIndex, integratedValue, ADD_VALUES); CHKERRV(ierr);
            
          }  // bComponent
        }  // bDof       
      }  // aComponent
    }  // aDof
  }  // elementNo
  
  // because this is used in nonlinear solver context, assembly has to be called here, not via data->finalAssembly
  MatAssemblyBegin(tangentStiffnessMatrix,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(tangentStiffnessMatrix,MAT_FINAL_ASSEMBLY);
  
  if (!tangentStiffnessMatrixInitialized_)
  {
    VLOG(3) << "tangent stiffness matrix before zero rows columns: " << PetscUtility::getStringMatrix(tangentStiffnessMatrix);
    VLOG(3) << "number dirichletIndices: " << dirichletIndices_.size();
    
    // zero rows and columns for which Dirichlet BC is set 
    MatZeroRowsColumns(tangentStiffnessMatrix, dirichletIndices_.size(), dirichletIndices_.data(), 1.0, PETSC_IGNORE, PETSC_IGNORE);
  
    VLOG(3) << "tangent stiffness matrix after zero rows columns: " << PetscUtility::getStringMatrix(tangentStiffnessMatrix);
    
    // set option that all insert/add operations to new nonzero locations will be discarded. This keeps the nonzero structure forever.
    // (The diagonal entries will be set to different values, but that doesn't matter because the Dirichlet values for updates are 0 and thus independent of the diagonal scaling (d*Δx=0 -> Δx=0 independent of d))
    MatSetOption(tangentStiffnessMatrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
    
    tangentStiffnessMatrixInitialized_ = true;
  }
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
Mat &FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
tangentStiffnessMatrix()
{
  return this->data_.tangentStiffnessMatrix();
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
Vec &FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
rightHandSide()
{
  return this->data_.rightHandSide().values();
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
Vec &FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
displacements()
{
  return this->data_.displacements().values();
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
setDisplacements(Vec &solverDisplacementVariable)
{
  // if the computation uses reduced vectors, expand to full vectors
  if (this->data_.computeWithReducedVectors())
  {
    this->expandVector(solverDisplacementVariable, this->data_.displacements().values());
  }
  else
  {
    VecCopy(solverDisplacementVariable, this->data_.displacements().values());
  }
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
computeInternalVirtualWork(Vec &resultVec)
{
  LOG(TRACE) << "computeInternalVirtualWork";
 
  // get pointer to mesh object
  std::shared_ptr<BasisOnMeshType> mesh = this->data_.mesh();
  
  const int D = BasisOnMeshType::dim();  // = 3
  //const int nDofs = mesh->nDofs();
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  const int nElements = mesh->nElements();
  const int nUnknowsPerElement = nDofsPerElement*D;    // 3 directions for displacements per dof
  
  // define shortcuts for quadrature
  typedef Quadrature::TensorProduct<D,QuadratureType> QuadratureDD;
  
  // define type to hold evaluations of integrand for stiffness matrix
  typedef std::array<double, nUnknowsPerElement> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureDD::numberEvaluations()
          > EvaluationsArrayType;     // evaluations[nGP^D](nUnknows,nUnknows)
  
  // setup arrays used for integration
  std::array<VecD<D>, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();
  EvaluationsArrayType evaluationsArray;
  
  PetscErrorCode ierr;
  
  // get memory from PETSc where to store result (not used)
  //const PetscScalar *result;
  //VecGetArray(resultVec, &result);
  
  // log file for PK2 stress
  std::ofstream pk2file("pk2.txt", std::ios::out | std::ios::binary);
  if (!pk2file.is_open())
    LOG(FATAL) << "could not open pk2.txt";
  pk2file << "#PK2 stress, xi" << std::endl;
  
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
      VecD<D> xi = samplingPoints[samplingPointIndex];
      
      // compute the 3xD jacobian of the parameter space to world space mapping
      Tensor2<D> jacobianMaterial = MathUtility::transformToDxD<D,D>(BasisOnMeshType::computeJacobian(geometryReferenceValues, xi));
      double jacobianDeterminant;
      Tensor2<D> inverseJacobianMaterial = MathUtility::computeInverse<D>(jacobianMaterial, jacobianDeterminant);
      // jacobianMaterial[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx
      // inverseJacobianMaterial[columnIdx][rowIdx] = dxi_rowIdx/dX_columnIdx because of inverse function theorem
      
      checkInverseIsCorrect<D>(jacobianMaterial, inverseJacobianMaterial, "jacobianMaterial");
      
      // F
      Tensor2<D> deformationGradient = this->computeDeformationGradient(displacementValues, inverseJacobianMaterial, xi);
      double deformationGradientDeterminant = MathUtility::computeDeterminant<D>(deformationGradient);  // J
      
      Tensor2<D> rightCauchyGreen = this->computeRightCauchyGreenTensor(deformationGradient);  // C = F^T*F
      
      checkSymmetry<D>(rightCauchyGreen, "C");
      
      double rightCauchyGreenDeterminant;   // J^2
      Tensor2<D> inverseRightCauchyGreen = MathUtility::computeSymmetricInverse<D>(rightCauchyGreen, rightCauchyGreenDeterminant);  // C^-1
      
      checkInverseIsCorrect<D>(rightCauchyGreen, inverseRightCauchyGreen, "rightCauchyGreen");
      
      // invariants
      std::array<double,3> invariants = this->computeInvariants(rightCauchyGreen, rightCauchyGreenDeterminant);  // I_1, I_2, I_3
      std::array<double,2> reducedInvariants = this->computeReducedInvariants(invariants, deformationGradientDeterminant); // Ibar_1, Ibar_2
  
      // artifical pressure p
      double artificialPressureTilde;
      const double artificialPressure = this->computeArtificialPressure(deformationGradientDeterminant, artificialPressureTilde);
      
      // Pk2 stress tensor S = S_vol + S_iso (p.234)
      Tensor2<D> fictitiousPK2Stress;
      Tensor2<D> pk2StressIsochoric;
      
      //! compute 2nd Piola-Kirchhoff stress tensor S = 2*dPsi/dC and the fictitious PK2 Stress Sbar
      Tensor2<D> PK2Stress = this->computePK2Stress(artificialPressure, rightCauchyGreen, inverseRightCauchyGreen, reducedInvariants, deformationGradientDeterminant, fictitiousPK2Stress, pk2StressIsochoric);      
      
      // check if PK2 is the same as from explicit formula for Mooney-Rivlin (p.249)
      if (VLOG_IS_ON(1))
      {
        this->checkFictitiousPK2Stress(fictitiousPK2Stress, rightCauchyGreen, deformationGradientDeterminant, reducedInvariants);
      }
      
      std::array<VecD<D>,nDofsPerElement> gradPhi = mesh->getGradPhi(xi);
      // (column-major storage) gradPhi[M][a] = dphi_M / dxi_a
      // gradPhi[column][row] = gradPhi[dofIndex][i] = dphi_dofIndex/dxi_i, columnIdx = dofIndex, rowIdx = which direction
      
      VLOG(2) << "";
      VLOG(2) << "element " << elementNo << " xi: " << xi;
      VLOG(2) << "  geometryReferenceValues: " << geometryReferenceValues;
      VLOG(2) << "  displacementValues: " << displacementValues;
      VLOG(2) << "  Jacobian: J_phi=" << jacobianMaterial;
      VLOG(2) << "  jacobianDeterminant: J=" << jacobianDeterminant;
      VLOG(2) << "  inverseJacobianMaterial: J_phi^-1=" << inverseJacobianMaterial;
      VLOG(2) << "  deformationGradient: F=" << deformationGradient;
      VLOG(2) << "  deformationGradientDeterminant: det F=" << deformationGradientDeterminant;
      VLOG(2) << "  rightCauchyGreen: C=" << rightCauchyGreen;
      VLOG(2) << "  rightCauchyGreenDeterminant: det C=" << rightCauchyGreenDeterminant;
      VLOG(2) << "  inverseRightCauchyGreen: C^-1=" << inverseRightCauchyGreen;
      VLOG(2) << "  invariants: I1,I2,I3: " << invariants;
      VLOG(2) << "  reducedInvariants: Ibar1, Ibar2: " << reducedInvariants;
      VLOG(2) << "  artificialPressure: p=" << artificialPressure << ", artificialPressureTilde: pTilde=" << artificialPressureTilde;
      VLOG(2) << "  fictitiousPK2Stress: Sbar=" << fictitiousPK2Stress;
      VLOG(2) << "  pk2StressIsochoric: S_iso=" << pk2StressIsochoric;
      VLOG(2) << "  PK2Stress: S=" << PK2Stress;
      VLOG(2) << "  gradPhi: " << gradPhi;
      
      VLOG(1) << "  xi: " << xi << ", J: " << deformationGradientDeterminant << ", p: " << artificialPressure << ", S11: " << PK2Stress[0][0];
      
      if (samplingPointIndex == 0 && D == 3)
        LOG(DEBUG) << " F11: " << deformationGradient[0][0] << ", F22,F33: " << deformationGradient[1][1] <<"," << deformationGradient[2][2]
          << ", F12,F13,F23: " << deformationGradient[1][0] << "," << deformationGradient[2][0] << "," << deformationGradient[2][1]
          << ", J: " << deformationGradientDeterminant << ", p: " << artificialPressure;
      
      pk2file << PK2Stress << " " << xi << std::endl;
      
      if (VLOG_IS_ON(2))
      {
        Tensor2<D> greenLangrangeStrain = this->computeGreenLagrangeStrain(rightCauchyGreen);
        VLOG(2) << "  strain E=" << greenLangrangeStrain;
      }
      
      // loop basis functions and evaluate integrand at xi
      for (int aDof = 0; aDof < nDofsPerElement; aDof++)           // index over dofs, each dof has D components, L in derivation
      {
        for (int aComponent = 0; aComponent < D; aComponent++)     // lower-case a in derivation, index over displacement components
        {
          // compute index of degree of freedom and component (result vector index)
          const int i = D*aDof + aComponent;
         
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
              double dphiL_dXA = 0.0;
              double dphiL_dXB = 0.0;
              
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
              
              integrand += 1./2. * PK2Stress[bInternal][aInternal] * (faB*dphiL_dXA + faA*dphiL_dXB);
                         
            }  // B, bInternal
          }  // A, aInternal
          
          // store integrand in evaluations array
          evaluationsArray[samplingPointIndex][i] = integrand;
          
        }  // a
      }  // L
    }  // function evaluations
    
    pk2file.close();
    
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
        const int i = D*aDof + aComponent;
   
        // integrate value and set entry in stiffness matrix
        double integratedValue = integratedValues[i];
        
        // Compute indices in stiffness matrix. For each dof there are 3 values for the 3 displacement components, even if D=2
        // Therefore the nDofsPerElement number is not the number of unknows.
        dof_no_t resultVectorIndex = dofNo[aDof]*3 + aComponent;
        
        //VLOG(2) << "  result vector "<<aDof<<","<<aComponent<<" = " <<i<<", dof " << dofNo[aDof];
        //VLOG(2) << "      vector index "<<resultVectorIndex<<", integrated value: "<<integratedValue;
            
        ierr = VecSetValue(resultVec, resultVectorIndex, integratedValue, INSERT_VALUES); CHKERRV(ierr);
        
      }  // aComponent
    }  // aDof
  }  // elementNo
  
  ierr = VecAssemblyBegin(resultVec); CHKERRV(ierr);
  ierr = VecAssemblyEnd(resultVec); CHKERRV(ierr);
  // return memory acces of result vector back to PETSc (not used)
  //VecRestoreArray(resultVec, &result);
  
  VLOG(1) << "    disp: " << PetscUtility::getStringVector(this->data_.displacements().values());
  VLOG(1) << "  ->wInt: " << PetscUtility::getStringVector(resultVec);
}


template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
computeExternalVirtualWork(Vec &resultVec)
{
  LOG(TRACE) << "computeExternalVirtualWork";
 
  // get pointer to mesh object
  std::shared_ptr<BasisOnMeshType> mesh = this->data_.mesh();
  
  const int D = BasisOnMeshType::dim();  // = 3
  //const int nDofs = mesh->nDofs();
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  //const int nElements = mesh->nElements();
  // define shortcuts for quadrature
  typedef Quadrature::TensorProduct<D,QuadratureType> QuadratureDD;
  typedef Quadrature::TensorProduct<D-1,QuadratureType> QuadratureSurface;
  
  // define type to hold evaluations of integrand for result vector, volume integrals for body forces
  typedef std::array<double, nDofsPerElement*D> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureDD::numberEvaluations()
          > EvaluationsArrayType;     // evaluations[nGP^D](nDofs*D)
  
  // setup arrays used for integration
  std::array<VecD<D>, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();
  EvaluationsArrayType evaluationsArray;
  
  // define type to hold evaluations of integrand for result vector, surface integrals for traction
  typedef std::array<
            EvaluationsType,
            QuadratureSurface::numberEvaluations()
          > EvaluationsArraySurfaceType;     // evaluations[nGP^D](nDofs*D)
  
  std::array<std::array<double,D-1>, QuadratureSurface::numberEvaluations()> samplingPointsSurface = QuadratureSurface::samplingPoints();
  EvaluationsArraySurfaceType evaluationsArraySurface;
  
  PetscErrorCode ierr;
  
  // zero result vector
  VecZeroEntries(resultVec);
  
 
  // -------------- body force in reference configuration --------------
  // loop over elements in bodyForceReferenceConfiguration_ that have a force vector specified
  for (typename std::vector<std::pair<element_no_t, VecD<D>>>::const_iterator iter = bodyForceReferenceConfiguration_.begin(); 
       iter != bodyForceReferenceConfiguration_.end(); iter++)
  {
    int elementNo = iter->first;
    VecD<D> forceVector = iter->second;
 
    // loop over integration points (e.g. gauss points)
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // get parameter values of current sampling point
      VecD<D> xi = samplingPoints[samplingPointIndex];
      
      // loop over dofs of element
      for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
      {
        VecD<D> dofIntegrand = forceVector * MathUtility::sqr(mesh->phi(dofIndex, xi));
        
        // store integrand in evaluations array
        for (int i = 0; i < D; i++)
        {
          evaluationsArray[samplingPointIndex][dofIndex*D+i] = dofIntegrand[i];
        }
        
      }  // dofIndex
          
    }  // function evaluations
    
    // integrate all values for result vector entries at once
    EvaluationsType integratedValues = QuadratureDD::computeIntegral(evaluationsArray);
    
    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNo = mesh->getElementDofNos(elementNo);
   
    // add entries in result vector
    // loop over indices of unknows (dofIndex,dofComponent)
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      for (int dofComponent = 0; dofComponent < D; dofComponent++)
      {
        // compute index of degree of freedom and component (integrade values vector index)
        const int i = D*dofIndex + dofComponent;
   
        // get integrated value
        double integratedValue = integratedValues[i];
        
        // Compute indices in result vector. For each dof there are 3 values for the 3 displacement components, even for D=2. 
        // Therefore the nDofsPerElement number is not the number of unknows.
        dof_no_t resultVectorIndex = dofNo[dofIndex]*3 + dofComponent;
        
        ierr = VecSetValue(resultVec, resultVectorIndex, integratedValue, ADD_VALUES); CHKERRV(ierr);
        
      }  // dofComponent
    }  // dofIndex
  }  // elementNo
  
  // -------------- body force in current configuration --------------
  // loop over elements in bodyForceCurrentConfiguration_ that have a force vector specified
  for (typename std::vector<std::pair<element_no_t, VecD<D>>>::const_iterator iter = bodyForceCurrentConfiguration_.begin(); 
       iter != bodyForceCurrentConfiguration_.end(); iter++)
  {
    int elementNo = iter->first;
    VecD<D> forceVector = iter->second;
    
    // get geometry field of reference configuration
    std::array<Vec3,nDofsPerElement> geometryReferenceValues;
    this->data_.geometryReference().getElementValues(elementNo, geometryReferenceValues);
    
    // loop over integration points (e.g. gauss points)
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // get parameter values of current sampling point
      VecD<D> xi = samplingPoints[samplingPointIndex];
      
      // compute the 3xD jacobian of the parameter space to world space mapping
      Tensor2<D> jacobianMaterial = MathUtility::transformToDxD<D,D>(BasisOnMeshType::computeJacobian(geometryReferenceValues, xi));
      double jacobianDeterminant = MathUtility::computeDeterminant<D>(jacobianMaterial);
      
      // loop over dofs of element
      for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
      {
        // scale the force vector with the Jacobian determinant: B = J*b
        // for incompressible material this should have no effect, J=1
 
        VecD<D> dofIntegrand = forceVector * jacobianDeterminant * MathUtility::sqr(mesh->phi(dofIndex, xi));
        
        // store integrand in evaluations array
        for (int i = 0; i < D; i++)
        {
          evaluationsArray[samplingPointIndex][dofIndex*D + i] = dofIntegrand[i];
        }
        
      }  // dofIndex
          
    }  // function evaluations
    
    // integrate all values for result vector entries at once
    EvaluationsType integratedValues = QuadratureDD::computeIntegral(evaluationsArray);
    
    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNo = mesh->getElementDofNos(elementNo);
   
    // add entries in result vector
    // loop over indices of unknows (dofIndex,dofComponent)
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      for (int dofComponent = 0; dofComponent < D; dofComponent++)
      {
        // compute index of degree of freedom and component (integrade values vector index)
        const int i = D*dofIndex + dofComponent;
   
        // get integrated value
        double integratedValue = integratedValues[i];
        
        // Compute indices in result vector. For each dof there are 3 values for the 3 displacement components, even for D=2.
        // Therefore the nDofsPerElement number is not the number of unknows.
        dof_no_t resultVectorIndex = dofNo[dofIndex]*3 + dofComponent;
        
        ierr = VecSetValue(resultVec, resultVectorIndex, integratedValue, ADD_VALUES); CHKERRV(ierr);
        
      }  // dofComponent
    }  // dofIndex
  }  // elementNo
  
  // -------------- surface traction in reference configuration --------------
  for (typename std::vector<TractionBoundaryCondition>::const_iterator iter = tractionReferenceConfiguration_.begin(); 
       iter != tractionReferenceConfiguration_.end(); iter++)
  {
    /* 
    element_no_t elementGlobalNo;
    
    Mesh::face_t face;
    std::vector<std::pair<dof_no_t, Vec3>> dofVectors;  //<element-local dof no, value>
    */
    element_no_t elementNo = iter->elementGlobalNo;
    LOG(DEBUG) << "tractionReferenceConfiguration_ on element " << elementNo;
    
    // loop over integration points (e.g. gauss points)
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPointsSurface.size(); samplingPointIndex++)
    {
      // get parameter values of current sampling point
      std::array<double,D-1> xiSurface = samplingPointsSurface[samplingPointIndex];
      VecD<D> xi = Mesh::getXiOnFace(iter->face, xiSurface);
      
      // set all entries to 0
      evaluationsArraySurface[samplingPointIndex] = {0.0};
      
      // compute the traction at xi by summing contributions from elemental dofs
      VecD<D> tractionReferenceConfiguration({0.0});
      for (typename std::vector<std::pair<dof_no_t, VecD<D>>>::const_iterator dofVectorsIter = iter->dofVectors.begin(); 
           dofVectorsIter != iter->dofVectors.end(); 
           dofVectorsIter++)
      {
        int dofIndex = dofVectorsIter->first;
        VecD<D> forceVector = dofVectorsIter->second;
        
        tractionReferenceConfiguration += forceVector * mesh->phi(dofIndex, xi);
      }
      
      // loop over dofs of element with given traction, those have potentially non-zero virtual external work contribution
      for (typename std::vector<std::pair<dof_no_t, VecD<D>>>::const_iterator dofVectorsIter = iter->dofVectors.begin();
           dofVectorsIter != iter->dofVectors.end();
           dofVectorsIter++)
      {
        int dofIndex = dofVectorsIter->first;
        VecD<D> forceVector = dofVectorsIter->second;
        
        VecD<D> dofIntegrand = tractionReferenceConfiguration * mesh->phi(dofIndex, xi);
        
        LOG(DEBUG) << "  dofIndex " << dofIndex << ", xi=" << xi << ", traction: " << tractionReferenceConfiguration[0] 
          << " phi = " << mesh->phi(dofIndex, xi);
      
        // store integrand in evaluations array
        for (int i = 0; i < D; i++)
        {
          evaluationsArraySurface[samplingPointIndex][dofIndex*D+i] = dofIntegrand[i];
        }
      }  // dofIndex
          
    }  // function evaluations
    
    // integrate all values for result vector entries at once
    EvaluationsType integratedValues = QuadratureSurface::computeIntegral(evaluationsArraySurface);
    
    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNo = mesh->getElementDofNos(elementNo);
   
    // add entries in result vector
    // loop over indices of unknows (dofIndex,dofComponent)
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      for (int dofComponent = 0; dofComponent < D; dofComponent++)
      {
        // compute index of degree of freedom and component (integrade values vector index)
        const int i = D*dofIndex + dofComponent;
   
        // get integrated value
        double integratedValue = integratedValues[i];
        
        LOG(DEBUG) << "  dof " << dofIndex << ", component " << dofComponent << " integrated value: " << integratedValue;
        
        // Compute indices in result vector. For each dof there are 3 values for the 3 displacement components, even for D=2.
        // Therefore the nDofsPerElement number is not the number of unknows.
        dof_no_t resultVectorIndex = dofNo[dofIndex]*3 + dofComponent;
        
        ierr = VecSetValue(resultVec, resultVectorIndex, integratedValue, ADD_VALUES); CHKERRV(ierr);
        
      }  // dofComponent
    }  // dofIndex
  }  // elementNo
  
  // -------------- surface traction in current configuration --------------
  for (typename std::vector<TractionBoundaryCondition>::const_iterator iter = tractionCurrentConfiguration_.begin(); 
       iter != tractionCurrentConfiguration_.end(); iter++)
  {
    /* 
    element_no_t elementGlobalNo;
    
    Mesh::face_t face;
    std::vector<std::pair<dof_no_t, Vec3>> dofVectors;  //<element-local dof no, value>
    */
    element_no_t elementNo = iter->elementGlobalNo;
    
    // get geometry field of reference configuration
    std::array<Vec3,nDofsPerElement> geometryReferenceValues;
    this->data_.geometryReference().getElementValues(elementNo, geometryReferenceValues);
    
    // get displacement field
    std::array<Vec3,nDofsPerElement> displacementValues;
    this->data_.displacements().getElementValues(elementNo, displacementValues);
    
    // loop over integration points (e.g. gauss points)
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPointsSurface.size(); samplingPointIndex++)
    {
      // get parameter values of current sampling point
      std::array<double,D-1> xiSurface = samplingPointsSurface[samplingPointIndex];
      VecD<D> xi = Mesh::getXiOnFace(iter->face, xiSurface);
      
      // compute the normal in reference configuration of the specified face at point xi
      VecD<D> normalReferenceConfiguration = MathUtility::transformToD<D,3>(mesh->getNormal(iter->face, geometryReferenceValues, xi));
        
      // compute the 3xD jacobian of the parameter space to world space mapping
      Tensor2<D> jacobianMaterial = MathUtility::transformToDxD<D,D>(BasisOnMeshType::computeJacobian(geometryReferenceValues, xi));
      double jacobianDeterminant;
      Tensor2<D> inverseJacobianMaterial = MathUtility::computeInverse<D>(jacobianMaterial, jacobianDeterminant);
      // jacobianMaterial[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx
      // inverseJacobianMaterial[columnIdx][rowIdx] = dxi_rowIdx/dX_columnIdx because of inverse function theorem
      
      // F
      Tensor2<D> deformationGradient = this->computeDeformationGradient(displacementValues, inverseJacobianMaterial, xi);
      
      // J*F^{-T}
      Tensor2<D> cofactorMatrix = MathUtility::computeCofactorMatrix<D>(deformationGradient);
      
      VecD<D> normalCurrentConfiguration = cofactorMatrix * normalReferenceConfiguration;
      double surfaceStretch = MathUtility::norm<D>(normalCurrentConfiguration);  // ds/dS, because |dS| = |normalReferenceConfiguration| = 1
      
      LOG(DEBUG) << "surfaceStretch: " << surfaceStretch;
      
      // T = t ds/dS (ds, dS are scalar)
      // Nansons formula: ds = J F^-T dS (ds, dS are vectors), ds/dS = |ds|/|dS|  = |J F^-T dS|/|dS|, construct |dS| = 1
       
      // set all entries to 0
      evaluationsArraySurface[samplingPointIndex] = {0.0};
      
      // compute the traction at xi by summing contributions from elemental dofs
      VecD<D> tractionCurrentConfiguration({0.0});
      for (typename std::vector<std::pair<dof_no_t, VecD<D>>>::const_iterator dofVectorsIter = iter->dofVectors.begin(); 
           dofVectorsIter != iter->dofVectors.end(); 
           dofVectorsIter++)
      {
        int dofIndex = dofVectorsIter->first;
        VecD<D> forceVector = dofVectorsIter->second;
        
        tractionCurrentConfiguration += forceVector * mesh->phi(dofIndex, xi);
      }
      
      // loop over dofs of element
      for (typename std::vector<std::pair<dof_no_t, VecD<D>>>::const_iterator dofVectorsIter = iter->dofVectors.begin();
           dofVectorsIter != iter->dofVectors.end(); 
           dofVectorsIter++)
      {
        int dofIndex = dofVectorsIter->first;
        
        VecD<D> dofIntegrand = tractionCurrentConfiguration * surfaceStretch * MathUtility::sqr(mesh->phi(dofIndex, xi));
        
        // store integrand in evaluations array
        for (int i = 0; i < D; i++)
        {
          evaluationsArraySurface[samplingPointIndex][dofIndex*D+i] = dofIntegrand[i];
        }
      }  // dofIndex
          
    }  // function evaluations
    
    // integrate all values for result vector entries at once
    EvaluationsType integratedValues = QuadratureSurface::computeIntegral(evaluationsArraySurface);
    
    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNo = mesh->getElementDofNos(elementNo);
   
    // add entries in result vector
    // loop over indices of unknows (dofIndex,dofComponent)
    for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
    {
      for (int dofComponent = 0; dofComponent < D; dofComponent++)
      {
        // compute index of degree of freedom and component (integrade values vector index)
        const int i = D*dofIndex + dofComponent;
   
        // get integrated value
        double integratedValue = integratedValues[i];
        
        // Compute indices in result vector. For each dof there are 3 values for the 3 displacement components, even for D=2.
        // Therefore the nDofsPerElement number is not the number of unknows.
        dof_no_t resultVectorIndex = dofNo[dofIndex]*3 + dofComponent;
        
        ierr = VecSetValue(resultVec, resultVectorIndex, integratedValue, ADD_VALUES); CHKERRV(ierr);
        
      }  // dofComponent
    }  // dofIndex
  }  // elementNo
  
  ierr = VecAssemblyBegin(resultVec); CHKERRV(ierr);
  ierr = VecAssemblyEnd(resultVec); CHKERRV(ierr);
  // return memory acces of result vector back to PETSc (not used)
  //VecRestoreArray(resultVec, &result);
  
  VLOG(1) << "  ->wExt: " << PetscUtility::getStringVector(resultVec);
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
getExternalVirtualWork(Vec &resultVec)
{
  if (this->data_.externalVirtualWorkIsConstant())
  {
    // the external virtual work is independent of the current displacements and was already computed from initializeBoundaryConditions
    VecCopy(this->data_.externalVirtualWork().values(), resultVec);
  }
  else
  {
    computeExternalVirtualWork(resultVec);
    this->data_.externalVirtualWork().values() = resultVec;
  }
  
}
 
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
computeInternalMinusExternalVirtualWork(Vec &resultVec)
{
  if (this->data_.computeWithReducedVectors())
  { 
    Vec &wInt = this->data_.internalVirtualWork().values();
    Vec &wExt = this->data_.externalVirtualWork().values();
    Vec &wIntReduced = this->data_.rhsReduced();
    
    computeInternalVirtualWork(wInt);
    getExternalVirtualWork(wExt);
    
    reduceVector(wInt, wIntReduced);
    reduceVector(wExt, resultVec);
   
    LOG(DEBUG) << "--                           dW_int: " << PetscUtility::getStringVector(wInt);
    LOG(DEBUG) << "--                           dW_ext: " << PetscUtility::getStringVector(wExt);
    
    // compute wInt - wExt
    PetscErrorCode ierr;
    ierr = VecAYPX(resultVec, -1, wIntReduced); CHKERRV(ierr);
  }
  else
  {    
    Vec &wInt = this->data_.internalVirtualWork().values();
    
    computeInternalVirtualWork(wInt);
    getExternalVirtualWork(resultVec);
    
    LOG(DEBUG) << "--                           dW_int: " << PetscUtility::getStringVector(wInt);
    LOG(DEBUG) << "--                           dW_ext: " << PetscUtility::getStringVector(resultVec);
    
    // compute wInt - wExt
    PetscErrorCode ierr;
    ierr = VecAYPX(resultVec, -1, wInt); CHKERRV(ierr);
  }
} 
 
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
initialize()
{
  initializeMaterialParameters();
  this->data_.initialize();
  
  // initialize external virtual energy
  initializeBoundaryConditions();
  printBoundaryConditions();
  
  //this->setStiffnessMatrix();
  //this->setRightHandSide();
  this->data_.finalAssembly();
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
initializeMaterialParameters()
{
  std::array<double,Term::nMaterialParameters> parameters 
    = PythonUtility::template getOptionArray<double, Term::nMaterialParameters>(this->specificSettings_, "materialParameters", 0.0);
 
  LOG(DEBUG) << "Material has " << Term::nMaterialParameters << " parameters, parsed from \"materialParameters\": " << parameters;
  
  std::vector<double> parametersVector(Term::nMaterialParameters);
  std::copy(parameters.begin(), parameters.end(), parametersVector.begin());
  
  // set all PARAM(i) values to the values given by materialParameters
  SEMT::set_parameters<Term::nMaterialParameters>::to(parametersVector);
}
 
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::TractionBoundaryCondition::
TractionBoundaryCondition(PyObject *specificSettings, std::shared_ptr<BasisOnMeshType> mesh)
{
  const int D = BasisOnMeshType::dim();
 
  // store global element no
  elementGlobalNo = PythonUtility::getOptionInt(specificSettings, "element", 0, PythonUtility::NonNegative);
  
  // store face
  std::string faceStr = PythonUtility::getOptionString(specificSettings, "face", "0+");
  face = Mesh::parseFace(faceStr);
  
  if (PythonUtility::hasKey(specificSettings, "constantValue") && PythonUtility::hasKey(specificSettings, "constantVector"))
  {
    LOG(ERROR) << "Specified both \"constantValue\" and \"constantVector\".";
  }
  
  // parse dof vectors
  if (PythonUtility::hasKey(specificSettings, "constantValue") || PythonUtility::hasKey(specificSettings, "constantVector"))
  {
    VecD<D> constantVector;
    
    if (PythonUtility::hasKey(specificSettings, "constantValue"))
    {
      double constantValue = PythonUtility::getOptionDouble(specificSettings, "constantValue", 0.0);
      
      VecD<D-1> xiSurface;
      for (int i = 0; i < D-1; i++)
        xiSurface[i] = 0.5;
       
      VecD<D> xi = Mesh::getXiOnFace(face, xiSurface);
      constantVector = MathUtility::transformToD<D,3>(mesh->getNormal(face, elementGlobalNo, xi) * constantValue);
    }
    else if (PythonUtility::hasKey(specificSettings, "constantVector"))
    {
      constantVector = PythonUtility::getOptionArray<double,D>(specificSettings, "constantVector", 0.0);
    }
    
    // get dofs indices within element that correspond to the selected face
    const int D = BasisOnMeshType::dim();
    const int nDofs = BasisOnMesh::BasisOnMeshBaseDim<D-1,typename BasisOnMeshType::BasisFunction>::nDofsPerElement();
    std::array<dof_no_t,nDofs> dofIndices;
    BasisOnMeshType::getFaceDofs(face, dofIndices);

    for (int i = 0; i < nDofs; i++)
    {
      dofVectors.push_back(std::pair<dof_no_t, VecD<D>>(dofIndices[i], constantVector));
    }
  }
  else if (PythonUtility::hasKey(specificSettings, "dofVectors"))
  {
    std::pair<dof_no_t, PyObject *> dofVectorItem;
    
    // loop over dofVectors
    for (dofVectorItem = PythonUtility::getOptionDictBegin<dof_no_t, PyObject *>(specificSettings, "dofVectors");
      !PythonUtility::getOptionDictEnd(specificSettings, "dofVectors"); 
      PythonUtility::getOptionDictNext<dof_no_t, PyObject *>(specificSettings, "dofVectors", dofVectorItem))
    {
      dof_no_t dofIndex = dofVectorItem.first;
      VecD<D> dofVector = PythonUtility::convertFromPython<double,D>(dofVectorItem.second);
      
      dofVectors.push_back(std::pair<dof_no_t, VecD<D>>(dofIndex, dofVector));
    }
  }
  else
  {
    LOG(ERROR) << "Traction on element has not specified \"constantValue\", \"constantVector\" nor \"dofVectors\".";
  }
}
  
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
initializeBoundaryConditions()
{
  LOG(TRACE)<<"initializeBoundaryConditions";
 
  dof_no_t nUnknowns = this->data_.nUnknowns();
  const int D = BasisOnMeshType::dim();
  Vec &externalVirtualWork = this->data_.externalVirtualWork().values();
  
  // ----- Neumann BC ----------
  // parse values for traction and body force
  this->data_.externalVirtualWorkIsConstant() = true;
  if (PythonUtility::hasKey(this->specificSettings_, "tractionReferenceConfiguration"))
  {
    PyObject *listItem = PythonUtility::getOptionListBegin<PyObject*>(this->specificSettings_, "tractionReferenceConfiguration"); 
    
    // loop over items in tractionReferenceConfiguration list
    for (;
         !PythonUtility::getOptionListEnd(this->specificSettings_, "tractionReferenceConfiguration"); 
         PythonUtility::getOptionListNext<PyObject*>(this->specificSettings_, "tractionReferenceConfiguration", listItem))
    {
      tractionReferenceConfiguration_.emplace_back(listItem, this->data_.mesh());
    }
  }
  
  if (PythonUtility::hasKey(this->specificSettings_, "tractionCurrentConfiguration"))
  {
    this->data_.externalVirtualWorkIsConstant() = false;
    
    PyObject *listItem = PythonUtility::getOptionListBegin<PyObject*>(this->specificSettings_, "tractionReferenceConfiguration"); 
    
    // loop over items in tractionCurrentConfiguration list
    for (; 
         !PythonUtility::getOptionListEnd(this->specificSettings_, "tractionCurrentConfiguration"); 
         PythonUtility::getOptionListNext<PyObject*>(this->specificSettings_, "tractionCurrentConfiguration", listItem))
    {
      tractionCurrentConfiguration_.emplace_back(listItem, this->data_.mesh());
    }
  }
  
  if (PythonUtility::hasKey(this->specificSettings_, "bodyForceReferenceConfiguration"))
  {
   // example entry: {0: [tmax,0,0], 5: [tmax,tmax,tmax]},   # {<element global no.>: <vector>, ...}
   // loop over dict items
   std::pair<element_no_t, PyObject *> dofVectorItem;
   for (dofVectorItem = PythonUtility::getOptionDictBegin<element_no_t, PyObject *>(this->specificSettings_, "bodyForceReferenceConfiguration");
      !PythonUtility::getOptionDictEnd(this->specificSettings_, "bodyForceReferenceConfiguration"); 
      PythonUtility::getOptionDictNext<element_no_t, PyObject *>(this->specificSettings_, "bodyForceReferenceConfiguration", dofVectorItem))
    {
      element_no_t elementGlobalNo = dofVectorItem.first;
      VecD<D> vector = PythonUtility::convertFromPython<VecD<D>>(dofVectorItem.second);
     
      bodyForceReferenceConfiguration_.push_back(std::pair<element_no_t,VecD<D>>(elementGlobalNo, vector));
    }
  }
  
  if (PythonUtility::hasKey(this->specificSettings_, "bodyForceCurrentConfiguration"))
  {
    this->data_.externalVirtualWorkIsConstant() = false;
     
    // example entry: {0: [tmax,0,0], 5: [tmax,tmax,tmax]},   # {<global dof no.>: <vector>, ...}
    // loop over dict items
    std::pair<element_no_t, PyObject *> dofVectorItem;
    for (dofVectorItem = PythonUtility::getOptionDictBegin<element_no_t, PyObject *>(this->specificSettings_, "bodyForceCurrentConfiguration");
       !PythonUtility::getOptionDictEnd(this->specificSettings_, "bodyForceCurrentConfiguration");
       PythonUtility::getOptionDictNext<element_no_t, PyObject *>(this->specificSettings_, "bodyForceCurrentConfiguration", dofVectorItem))
    {
      element_no_t elementGlobalNo = dofVectorItem.first;
      VecD<D> vector = PythonUtility::convertFromPython<VecD<D>>(dofVectorItem.second);
     
      bodyForceReferenceConfiguration_.push_back(std::pair<element_no_t,VecD<D>>(elementGlobalNo, vector));
    }
  }
  
  // if the external virtual work is constant, i.e. independent of the displacement, it does not have to be computed for every new iteration of the nonlinear solver. Compute it once, now.
  if (this->data_.externalVirtualWorkIsConstant())
  {
    computeExternalVirtualWork(externalVirtualWork);
    
    VLOG(1) << "-- initially computed dW_ext: " << PetscUtility::getStringVector(externalVirtualWork);
    VLOG(1) << "-- initially computed dW_ext: " << PetscUtility::getStringVector(this->data_.externalVirtualWork().values());
  }
  
  // --------- Dirichlet BC --------------
  // iterate over Dirichlet boundary conditions
  
  // get the first dirichlet boundary condition from the list
  std::pair<node_no_t, double> boundaryCondition 
    = PythonUtility::getOptionDictBegin<dof_no_t, double>(this->specificSettings_, "dirichletBoundaryCondition");
  
  // loop over Dirichlet boundary conditions
  for (; !PythonUtility::getOptionDictEnd(this->specificSettings_, "dirichletBoundaryCondition"); 
       PythonUtility::getOptionDictNext<dof_no_t, double>(this->specificSettings_, "dirichletBoundaryCondition", boundaryCondition))
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
    
    // store values
    dirichletIndices_.push_back(boundaryConditionUnknownsIndex);
    dirichletValues_.push_back(boundaryConditionValue);
  }
  zeros_.resize(dirichletValues_.size(), 0.0);
  
  const int vectorSize = this->data_.mesh()->nDofs() * D - dirichletValues_.size();   // dofs always contain 3 entries for every entry (x,y,z)
  this->data_.initializeReducedVariables(vectorSize);
  
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
printBoundaryConditions()
{
  const int D = BasisOnMeshType::dim();
  /*
   * 
  std::vector<dof_no_t> dirichletIndices_;  ///< the indices of unknowns (not dofs) for which the displacement is fixed
  std::vector<double> dirichletValues_;     ///< the to dirichletIndices corresponding fixed values for the displacement
  std::vector<double> rhsValues_;           ///< the entries in rhs for which for u Dirichlet values are specified
  
  //TODO split into boundary conditions class
  struct TractionBoundaryCondition
  {
    element_no_t elementGlobalNo;
    
    Mesh::face_t face;
    std::vector<std::pair<dof_no_t, Vec3>> dofVectors;  //<element-local dof no, value>
    
    // parse values from python config, e.g. {"element": 1, "face": "0+", "dofVectors:", {0: [tmax,0,0], 1: [tmax,0,0], 2: [tmax,0,0], 3: [tmax,0,0]}}
    TractionBoundaryCondition(PyObject *specificSettings);
  };
  
  std::vector<TractionBoundaryCondition> tractionReferenceConfiguration_;   //< tractions for elements
  std::vector<TractionBoundaryCondition> tractionCurrentConfiguration_;
  
  std::vector<std::pair<dof_no_t, Vec3>> bodyForceReferenceConfiguration_;  //< <global dof no, vector>
  std::vector<std::pair<dof_no_t, Vec3>> bodyForceCurrentConfiguration_;    //< <global dof no, vector>
  */
  LOG(DEBUG) << "============ Dirichlet BC ============";
  LOG(DEBUG) << "dirichletIndices_: " << dirichletIndices_;
  LOG(DEBUG) << "dirichletValues_: " << dirichletValues_ << std::endl;
  
  LOG(DEBUG) << "============ Neumann BC ============";
  LOG(DEBUG) << "tractionReferenceConfiguration_: size: " << tractionReferenceConfiguration_.size();
  for (int i = 0; i < tractionReferenceConfiguration_.size(); i++)
  {
    LOG(DEBUG) << "no. " << i << ", element " << tractionReferenceConfiguration_[i].elementGlobalNo 
      << ", face: " << tractionReferenceConfiguration_[i].face << ", " << tractionReferenceConfiguration_[i].dofVectors.size() << " dofVectors: ";
    for (int j = 0; j < tractionReferenceConfiguration_[i].dofVectors.size(); j++)
    {
      std::pair<dof_no_t, VecD<D>> &dofVector = tractionReferenceConfiguration_[i].dofVectors[j];
      LOG(DEBUG) << "   " << j << ". dof index: " << dofVector.first << ", vector: " << dofVector.second;
    }
  }
  
  LOG(DEBUG) << "tractionCurrentConfiguration_: size: " << tractionCurrentConfiguration_.size();
  for (int i = 0; i < tractionCurrentConfiguration_.size(); i++)
  {
    LOG(DEBUG) << "no. " << i << ", element " << tractionCurrentConfiguration_[i].elementGlobalNo 
      << ", face: " << tractionCurrentConfiguration_[i].face << ", " << tractionCurrentConfiguration_[i].dofVectors.size() << " dofVectors: ";
    for (int j = 0; j < tractionCurrentConfiguration_[i].dofVectors.size(); j++)
    {
      std::pair<dof_no_t, VecD<D>> &dofVector = tractionCurrentConfiguration_[i].dofVectors[j];
      LOG(DEBUG) << "   " << j << ". dof index: " << dofVector.first << ", vector: " << dofVector.second;
    }
  }
  LOG(DEBUG) << "bodyForceReferenceConfiguration_: size: " << bodyForceReferenceConfiguration_.size();
  
  for (int i = 0; i < bodyForceReferenceConfiguration_.size(); i++)
  {
    LOG(DEBUG) << "   " << i << ". dof global no: " << bodyForceReferenceConfiguration_[i].first
      << ", vector: " << bodyForceReferenceConfiguration_[i].second;
  }
      
  LOG(DEBUG) << "bodyForceCurrentConfiguration_: size: " << bodyForceCurrentConfiguration_.size();
  
  for (int i = 0; i < bodyForceCurrentConfiguration_.size(); i++)
  {
    LOG(DEBUG) << "   " << i << ". dof global no: " << bodyForceCurrentConfiguration_[i].first
      << ", vector: " << bodyForceCurrentConfiguration_[i].second;
  }
  LOG(DEBUG) << "============ ============";
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
reduceVector(Vec &input, Vec &output)
{
  // compute number of reduced dofs
  const int nDofs = this->data_.mesh()->nDofs() * 3;
  int nDofsReduced = nDofs - this->dirichletValues_.size();
  
  // for 2D problems fix 3rd displacement values to 0
  if (BasisOnMeshType::dim() == 2)
  {
    nDofsReduced -= this->data_.mesh()->nDofs();
  }
  
  const double *inputData;
  double *outputData;
  VecGetArrayRead(input, &inputData);
  VecGetArray(output, &outputData);
  
  dof_no_t reducedIndex = 0;
  std::vector<dof_no_t>::const_iterator dirichletIndicesIter = dirichletIndices_.begin();
  
  for (dof_no_t currentDofNo = 0; currentDofNo < nDofs; currentDofNo++)
  {
    // exclude variables for which Dirichlet BC are set
    if (currentDofNo == *dirichletIndicesIter)
    {
      dirichletIndicesIter++;
      continue;
    }
    
    // for 2D problems exclude every 3rd entry that corresponds to z-displacements
    if (BasisOnMeshType::dim() == 2 && currentDofNo % 3 == 2)
    {
      continue;
    }
    
    outputData[reducedIndex++] = inputData[currentDofNo];
  }
  
  VecRestoreArrayRead(input, &inputData);
  VecRestoreArray(output, &outputData);
}


template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
expandVector(Vec &input, Vec &output)
{
  // compute number of reduced dofs
  const int nDofs = this->data_.mesh()->nDofs() * 3;
  int nDofsReduced = nDofs - this->dirichletValues_.size();
  
  // for 2D problems fix 3rd displacement values to 0
  if (BasisOnMeshType::dim() == 2)
  {
    nDofsReduced -= this->data_.mesh()->nDofs();
  }
  
  const double *inputData;
  double *outputData;
  VecGetArrayRead(input, &inputData);
  VecGetArray(output, &outputData);
  
  dof_no_t reducedIndex = 0;
  std::vector<dof_no_t>::const_iterator dirichletIndicesIter = dirichletIndices_.begin();
  std::vector<double>::const_iterator dirichletValuesIter = dirichletValues_.begin();
  
  for (dof_no_t currentDofNo = 0; currentDofNo < nDofs; currentDofNo++)
  {
    // exclude variables for which Dirichlet BC are set
    if (currentDofNo == *dirichletIndicesIter)
    {
      outputData[currentDofNo] = *dirichletValuesIter;
      dirichletIndicesIter++;
      dirichletValuesIter++;
      
      continue;
    }
    
    // for 2D problems exclude every 3rd entry that corresponds to z-displacements
    if (BasisOnMeshType::dim() == 2 && currentDofNo % 3 == 2)
    {
      outputData[currentDofNo] = 0.0;
      continue;
    }
    
    outputData[currentDofNo] = inputData[reducedIndex++];
  }
  
  VecRestoreArrayRead(input, &inputData);
  VecRestoreArray(output, &outputData);
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
applyDirichletBoundaryConditionsInNonlinearFunction(Vec &f)
{
  LOG(TRACE) << "applyDirichletBoundaryConditionsInNonlinearFunction";
  
  if (!this->data_.computeWithReducedVectors())
  {
    // overwrite the values in f that correspond to entries for which Dirichlet BC are set with zero values,
    // such that the nonlinear equation is satisfied for them
    PetscErrorCode ierr;
    ierr = VecSetValues(f, dirichletIndices_.size(), dirichletIndices_.data(), zeros_.data(), INSERT_VALUES); CHKERRV(ierr);
  }
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
applyDirichletBoundaryConditionsInDisplacements()
{
  LOG(TRACE) << "applyDirichletBoundaryConditionsInDisplacements";
 
  PetscErrorCode ierr;
  
  // if the nonlinear solver should work with not reduced vectors
  if (!this->data_.computeWithReducedVectors())
  {
    // set entries of Dirichlet BCs to specified values
    ierr = VecSetValues(this->data_.displacements().values(), dirichletIndices_.size(), dirichletIndices_.data(), dirichletValues_.data(), INSERT_VALUES); CHKERRV(ierr);
    
    // for 2D problems fix 3rd displacement values to 0
    if (BasisOnMeshType::dim() == 2)
    {
      const int nDofs = this->data_.mesh()->nDofs() * 3;
      std::vector<int> indices;
      std::vector<double> zeros(nDofs, 0.0);
      
      indices.reserve(nDofs);
      for (int i = 2; i < nDofs; i += 3)
      {
        indices.push_back(i);
      }
      VLOG(1) << "(indices: " << indices << ")";
      ierr = VecSetValues(this->data_.displacements().values(), indices.size(), indices.data(), zeros.data(), INSERT_VALUES); CHKERRV(ierr);  
    }
    
    VLOG(1) << "after applying Dirichlet BC displacements u:" << PetscUtility::getStringVector(this->data_.displacements().values());
  }
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
applyDirichletBoundaryConditionsInStiffnessMatrix(Mat &matrix)
{
  // if the nonlinear solver should work with not reduced vectors
  if (!this->data_.computeWithReducedVectors())
  {
    // zero rows and columns for which Dirichlet BC is set 
    PetscErrorCode ierr;
    ierr = MatZeroRowsColumns(matrix, dirichletIndices_.size(), dirichletIndices_.data(), 1.0, PETSC_IGNORE, PETSC_IGNORE); CHKERRV(ierr);
  
    // for 2D problems clear also every 3rd entry which corresponds to z displacements
    if (BasisOnMeshType::dim() == 2)
    {
      std::vector<int> indices(this->data_.mesh()->nDofs());
      for (dof_no_t i = 0; i < this->data_.mesh()->nDofs(); i++)
      {
        indices[i] = i*3 + 2;
      }
      ierr = MatZeroRowsColumns(matrix, this->data_.mesh()->nDofs(), indices.data(), 1.0, PETSC_IGNORE, PETSC_IGNORE); CHKERRV(ierr);
    }
  }
}
  
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
reduceMatrix(Mat &input, Mat &output)
{
  // TODO
  //std::vector<int> 
 
  //ierr = ISCreateGeneral(PETSC_COMM_WORLD,PetscInt n,const PetscInt idx[],PetscCopyMode mode,IS *is)
  //ierr = MatGetSubMatrix(input, isrow, iscol, MAT_REUSE_MATRIX, &output); CHKERRV(ierr);
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
computeAnalyticalStiffnessMatrix(Mat &solverStiffnessMatrix)
{
  if (this->data_.computeWithReducedVectors())
  {
    // Compute new stiffness matrix from current displacements. This uses the full vectors.
    this->setStiffnessMatrix();
    
    // create new reduced stiffness matrix by copying the non-zero entries
    reduceMatrix(this->data_.tangentStiffnessMatrix(), solverStiffnessMatrix);
  }
  else 
  {
    // Compute new stiffness matrix from current displacements. This uses the full vectors.
    this->setStiffnessMatrix(solverStiffnessMatrix);
    
    // set the appropriate entries to match Dirichlet BC
    applyDirichletBoundaryConditionsInStiffnessMatrix(solverStiffnessMatrix);
  }
}
  
};    // namespace