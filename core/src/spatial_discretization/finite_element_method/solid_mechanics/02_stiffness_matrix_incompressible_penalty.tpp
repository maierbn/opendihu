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
  std::array<std::array<double,D>, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();
  EvaluationsArrayType evaluationsArray;
  
  typedef std::array<Vec3,BasisOnMeshType::dim()> Tensor2;
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
      ElasticityTensor elasticity = this->computeElasticityTensor(artificialPressure, artificialPressureTilde, rightCauchyGreen, inverseRightCauchyGreen, fictitiousPK2Stress, pk2StressIsochoric, deformationGradientDeterminant, reducedInvariants);

      std::array<Vec3,nDofsPerElement> gradPhi = mesh->getGradPhi(xi);
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
setDisplacements(Vec &displacements)
{
  VecCopy(displacements, this->data_.displacements().values());
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
  std::array<std::array<double,D>, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();
  EvaluationsArrayType evaluationsArray;
  
  typedef std::array<Vec3,BasisOnMeshType::dim()> Tensor2;
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
      std::array<double,D> xi = samplingPoints[samplingPointIndex];
      
      // compute the 3xD jacobian of the parameter space to world space mapping
      Tensor2 jacobianMaterial = BasisOnMeshType::computeJacobian(geometryReferenceValues, xi);
      double jacobianDeterminant;
      Tensor2 inverseJacobianMaterial = MathUtility::computeInverse(jacobianMaterial, jacobianDeterminant);
      // jacobianMaterial[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx
      // inverseJacobianMaterial[columnIdx][rowIdx] = dxi_rowIdx/dX_columnIdx because of inverse function theorem
      
      checkInverseIsCorrect(jacobianMaterial, inverseJacobianMaterial, "jacobianMaterial");
      
      // F
      Tensor2 deformationGradient = this->computeDeformationGradient(displacementValues, inverseJacobianMaterial, xi);
      double deformationGradientDeterminant = MathUtility::computeDeterminant(deformationGradient);  // J
      
      Tensor2 rightCauchyGreen = this->computeRightCauchyGreenTensor(deformationGradient);  // C = F^T*F
      
      checkSymmetry(rightCauchyGreen, "C");
      
      double rightCauchyGreenDeterminant;   // J^2
      Tensor2 inverseRightCauchyGreen = MathUtility::computeSymmetricInverse(rightCauchyGreen, rightCauchyGreenDeterminant);  // C^-1
      
      checkInverseIsCorrect(rightCauchyGreen, inverseRightCauchyGreen, "rightCauchyGreen");
      
      // invariants
      std::array<double,3> invariants = this->computeInvariants(rightCauchyGreen, rightCauchyGreenDeterminant);  // I_1, I_2, I_3
      std::array<double,2> reducedInvariants = this->computeReducedInvariants(invariants, deformationGradientDeterminant); // Ibar_1, Ibar_2
  
      // artifical pressure p
      double artificialPressureTilde;
      const double artificialPressure = this->computeArtificialPressure(deformationGradientDeterminant, artificialPressureTilde);
      
      // Pk2 stress tensor S = S_vol + S_iso (p.234)
      std::array<Vec3,3> fictitiousPK2Stress;
      std::array<Vec3,3> pk2StressIsochoric;
      
      //! compute 2nd Piola-Kirchhoff stress tensor S = 2*dPsi/dC and the fictitious PK2 Stress Sbar
      Tensor2 PK2Stress = this->computePK2Stress(artificialPressure, rightCauchyGreen, inverseRightCauchyGreen, reducedInvariants, deformationGradientDeterminant, fictitiousPK2Stress, pk2StressIsochoric);      
      
      // check if PK2 is the same as from explicit formula for Mooney-Rivlin (p.249)
      if (VLOG_IS_ON(1))
      {
        this->checkFictitiousPK2Stress(fictitiousPK2Stress, rightCauchyGreen, deformationGradientDeterminant, reducedInvariants);
      }
      
      std::array<Vec3,nDofsPerElement> gradPhi = mesh->getGradPhi(xi);
      // (column-major storage) gradPhi[M][a] = dphi_M / dxi_a
      // gradPhi[column][row] = gradPhi[dofIndex][i] = dphi_dofIndex/dxi_i, columnIdx = dofIndex, rowIdx = which direction
      
      VLOG(2) << "";
      VLOG(2) << "element " << elementNo << " xi: " << xi;
      VLOG(2) << "  geometryReferenceValues: " << geometryReferenceValues;
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
      
      pk2file << PK2Stress << " " << xi << std::endl;
      
      if (VLOG_IS_ON(2))
      {
        Tensor2 greenLangrangeStrain = this->computeGreenLagrangeStrain(rightCauchyGreen);
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
              
              integrand += PK2Stress[bInternal][aInternal] * (faB*dphiL_dXA + faA*dphiL_dXB);
                         
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
        
        // Compute indices in stiffness matrix. For each dof there are D=3 values for the 3 displacement components. 
        // Therefore the nDofsPerElement number is not the number of unknows.
        dof_no_t resultVectorIndex = dofNo[aDof]*D + aComponent;
        
        //VLOG(2) << "  result vector "<<aDof<<","<<aComponent<<" = " <<i<<", dof " << dofNo[aDof];
        //VLOG(2) << "      vector index "<<resultVectorIndex<<", integrated value: "<<integratedValue;
            
        ierr = VecSetValue(resultVec, resultVectorIndex, integratedValue, INSERT_VALUES); CHKERRV(ierr);
        
      }  // aComponent
    }  // aDof
  }  // elementNo
  
  VecAssemblyBegin(resultVec);
  VecAssemblyEnd(resultVec);
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
initialize()
{
  initializeMaterialParameters();
  this->data_.initialize();
  initializeBoundaryConditions();
  //this->setStiffnessMatrix();
  this->setRightHandSide();
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
  Vec &rightHandSide = this->data_.externalVirtualEnergy().values();
  
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
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
applyDirichletBoundaryConditionsInNonlinearFunction(Vec &f)
{
  LOG(TRACE) << "applyDirichletBoundaryConditionsInNonlinearFunction";

  // overwrite the values in f that correspond to entries for which Dirichlet BC are set with rhs values, such that nonlinear equation is satisfied for them
  VecSetValues(f, dirichletIndices_.size(), dirichletIndices_.data(), rhsValues_.data(), INSERT_VALUES);
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
 
  // set entries of Dirichlet BCs to specified values
  VecSetValues(this->data_.displacements().values(), dirichletIndices_.size(), dirichletIndices_.data(), dirichletValues_.data(), INSERT_VALUES);
  
  VLOG(1) << "after applying Dirichlet BC, displacements u:" << this->data_.displacements().values();
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
  // zero rows and columns for which Dirichlet BC is set 
  MatZeroRowsColumns(matrix, dirichletIndices_.size(), dirichletIndices_.data(), 1.0, PETSC_IGNORE, PETSC_IGNORE);
}
  
};    // namespace