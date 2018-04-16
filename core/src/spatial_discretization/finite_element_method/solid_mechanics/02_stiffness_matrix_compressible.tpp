#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_compressible.h"

#include <Python.h>  // has to be the first included header
#include <memory>
#include <petscksp.h>
#include <petscsys.h>

#include "semt/Semt.h"
#include "semt/Shortcuts.h"
#include "control/types.h"
#include "utility/math_utility.h"

namespace SpatialDiscretization
{
 
/* 
// initialization for mooney rivlin
template<typename BasisOnMeshType, typename MixedQuadratureType>
void FiniteElementMethodStiffnessMatrix<BasisOnMeshType, MixedQuadratureType, Equation::Static::CompressibleMooneyRivlin, Mesh::isDeformable<typename BasisOnMeshType::Mesh>>:: 
initialize()
{
  kappa_ = PythonUtility::getOptionDouble(this->specificSettings_, "kappa", 1000, PythonUtility::Positive);
  std::array<double,2> materialConstants = PythonUtility::getOptionArray<double,2>(this->specificSettings_, "c");
  
  c1_ = materialConstants[0];
  c2_ = materialConstants[1];
  FiniteElementMethodBase<BasisOnMeshType,MixedQuadratureType>::initialize();
}

// stiffness matrix for compressible mooney rivlin
template<typename BasisOnMeshType, typename MixedQuadratureType>
void FiniteElementMethodStiffnessMatrix<BasisOnMeshType, MixedQuadratureType, Equation::Static::CompressibleMooneyRivlin, Mesh::isDeformable<typename BasisOnMeshType::Mesh>>:: 
setStiffnessMatrix()
{
  LOG(TRACE)<<"setStiffnessMatrix for solid mechanics";

  // naming: P = pressure variables, U = displacement variables
  
  typedef typename BasisOnMeshType::HighOrderBasisOnMesh HighOrderBasisOnMesh;
  typedef typename BasisOnMeshType::LowOrderBasisOnMesh LowOrderBasisOnMesh;
  
  // get references to mesh objects
  std::shared_ptr<HighOrderBasisOnMesh> basisOnMeshU = this->data_.mixedMesh()->highOrderBasisOnMesh();
  std::shared_ptr<LowOrderBasisOnMesh> basisOnMeshP = this->data_.mixedMesh()->lowOrderBasisOnMesh();
  
  const int D = BasisOnMeshType::dim();
  //const int nDofsU = basisOnMeshU->nDofs();
  //const int nDofsP = basisOnMeshP->nDofs();
  const int nDofsUPerElement = HighOrderBasisOnMesh::nDofsPerElement();
  const int nDofsPPerElement = LowOrderBasisOnMesh::nDofsPerElement();
  const int nElements = basisOnMeshU->nElements();
  
  // define shortcuts for integrator and basis
  typedef Quadrature::TensorProduct<D,typename MixedQuadratureType::HighOrderQuadrature> QuadratureU;
  typedef Quadrature::TensorProduct<D,typename MixedQuadratureType::LowOrderQuadrature> QuadratureP;
  
  typedef MathUtility::Matrix<nDofsUPerElement,nDofsUPerElement> EvaluationsUType;
  typedef MathUtility::Matrix<nDofsUPerElement,nDofsPPerElement> EvaluationsUPType;
  typedef MathUtility::Matrix<nDofsPPerElement,nDofsPPerElement> EvaluationsPType;
  typedef std::array<
            EvaluationsUType,
            QuadratureU::numberEvaluations()
          > EvaluationsUUMatrixType;     // evaluations[nGP^D][nDofsUPerElement][nDofsUPerElement]
  typedef std::array<
            EvaluationsUPType,
            QuadratureU::numberEvaluations()
          > EvaluationsUPMatrixType;     // evaluations[nGP^D][nDofsPPerElement][nDofsPPerElement]
  typedef std::array<
            EvaluationsPType,
            QuadratureP::numberEvaluations()
          > EvaluationsPPMatrixType;     // evaluations[nGP^D][nDofsUPerElement][nDofsUPerElement]
  
  // setup arrays used for integration
  std::array<std::array<double,D>, QuadratureU::numberEvaluations()> samplingPointsU = QuadratureU::samplingPoints();
  std::array<std::array<double,D>, QuadratureP::numberEvaluations()> samplingPointsP = QuadratureP::samplingPoints();
  EvaluationsUUMatrixType evaluationsKUU({0});
  EvaluationsUPMatrixType evaluationsKUP({0});
  EvaluationsPPMatrixType evaluationsKPP({0});
  
  std::array<std::array<double,nDofsUPerElement>,QuadratureU::numberEvaluations()> evaluationsFU({0});
  std::array<std::array<double,nDofsPPerElement>,QuadratureP::numberEvaluations()> evaluationsFP({0});
  
  typedef std::array<Vec3,BasisOnMeshType::dim()> Tensor2;
  typedef std::array<double,21> Tensor4;   // data type for 4th order elasticity tensor with due to symmetry has 21 independent components
  
  // the matrix formulation follows the paper "A finitie element formultaion for nonlinear incomporessible elastic and inelastic analysis" by Sussmann and Bathe
  
  
  PetscErrorCode ierr;
  // loop over elements 
  for (int elementNo = 0; elementNo < nElements; elementNo++)
  {
    // get indices of element-local dofs
    std::vector<dof_no_t> dofUNo;
    basisOnMeshU->getElementDofNos(elementNo, dofUNo);
    //std::vector<dof_no_t> dofPNo = basisOnMeshP->getElementDofNos(elementNo);
    
    // get geometry field of current configuration of meshU
    std::array<Vec3,HighOrderBasisOnMesh::nDofsPerElement()> geometryCurrentElementalDofValues;
    basisOnMeshU->getElementGeometry(elementNo, geometryCurrentElementalDofValues);
        
    // get geometry field of reference configuration
    std::array<Vec3,HighOrderBasisOnMesh::nDofsPerElement()> geometryReferenceElementalDofValues;
    this->data_.geometryReference().template getElementValues<D>(elementNo, geometryReferenceElementalDofValues);
    
    // get displacement field
    std::array<Vec3,HighOrderBasisOnMesh::nDofsPerElement()> displacementElementalDofValues;
    this->data_.displacement().template getElementValues<D>(elementNo, displacementElementalDofValues);
    
    // get separately interpolated pressure field
    std::array<double,LowOrderBasisOnMesh::nDofsPerElement()> separatePressureElementalDofValues;
    this->data_.pressure().template getElementValues<1>(elementNo, separatePressureElementalDofValues);
    
    // get right hand side component for displacements in weak form, fu
    std::array<double,HighOrderBasisOnMesh::nDofsPerElement()> fValues;
    this->data_.f().template getElementValues<1>(elementNo, fValues);
    
    // loop over integration points (e.g. gauss points) for displacement field
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPointsU.size(); samplingPointIndex++)
    {
      // get parameter values of current sampling point
      std::array<double,D> xi = samplingPointsU[samplingPointIndex];
      
      // compute the 3xD jacobian of the parameter space to world space mapping
      Tensor2 jacobian = HighOrderBasisOnMesh::computeJacobian(geometryCurrent, xi);
      
      // F
      Tensor2 deformationGradient = this->computeDeformationGradient(geometryReference, displacement, jacobian, xi);
      double deformationGradientDeterminant = MathUtility::computeDeterminant(deformationGradient);  // J
      
      Tensor2 rightCauchyGreen = this->computeRightCauchyGreenTensor(deformationGradient);  // C = F^T*F
      
      double rightCauchyGreenDeterminant;   // J^2
      Tensor2 inverseRightCauchyGreen = MathUtility::computeSymmetricInverse(rightCauchyGreen, rightCauchyGreenDeterminant);  // C^-1
      
      std::array<double,3> invariants = this->computeInvariants(rightCauchyGreen, rightCauchyGreenDeterminant);  // I_1, I_2, I_3
      std::array<double,3> reducedInvariants = this->computeReducedInvariants(invariants, deformationGradientDeterminant);  // J_1, J_2, J_3
      const double J3 = reducedInvariants[2];
      
      // Pk2 stress tensor S = 
      Tensor2 PK2Stress = computePK2Stress(rightCauchyGreen, inverseRightCauchyGreen, invariants, reducedInvariants);
      
      // elasticity tensor C_{ijkl}
      Tensor4 elasticity = computeElasticityTensor(rightCauchyGreen, inverseRightCauchyGreen, invariants, reducedInvariants);
      
      // get the pressure from the separate field variable
      std::array<double,1> interpolatedPressure = LowOrderBasisOnMesh::interpolateValueInElement(separatePressureElementalDofValues, xi);
      
      double pressureFromField = interpolatedPressure[0];
      double pressureFromDisplacements = -kappa_*(J3 - 1);
      double pressureDifference = pressureFromDisplacements - pressureFromField;
                  
      double pp = kappa_;
      double wh = 1./2 * kappa_ * MathUtility::sqr(J3 - 1);
      
      double cpp = -1./pp;
      
      // compute evaluation at current quadrature sampling point for element stiffness matrices kuu, kup and kpp
      // row index
      for (int i = 0; i < nDofsUPerElement; i++)
      {
        int dofIndexI = i;
        
        // dp_tilde/dp_i where the first p_tilde is the interpolated pressure and p_i is the pressure DOF
        double dpdpi = LowOrderBasisOnMesh::phi(dofIndexI, xi);
          
        // column index 
        for (int j = 0; j < nDofsUPerElement; j++)
        {
          int dofIndexJ = j;
          
          // dp_tilde/dp_j where the first p_tilde is the interpolated pressure and p_j is the pressure DOF
          double dpdpj = LowOrderBasisOnMesh::phi(dofIndexJ, xi);
          
          // indices of elasticity tensor
          for (int k = 0; k < 3; k++)
          {
            for (int l = 0; l < 3; l++)
            {  
              // dp/deps_kl: derivative of \bar{p} w.r.t Green-Lagrange strain tensor E (or eps in paper)
              double dp_deps_kl = -kappa_*J3*rightCauchyGreen[l][k];
             
              // deps_kl/du_i
              double depskl_dui = 1./2*(rightCauchyGreen[l][i]*HighOrderBasisOnMesh::dphi_dxi(dofIndexI, k, xi)
                  + rightCauchyGreen[k][i]*HighOrderBasisOnMesh::dphi_dxi(dofIndexI, l, xi));
              
              // d^2eps_kl/(dui*duj)
              double d2epskl_duij = 0;
              if (i == j)   // d^2eps_kl/(dui*duj) is zero for i!=j
              {
                d2epskl_duij = 1./2*(HighOrderBasisOnMesh::dphi_dxi(dofIndexI, k, xi)*HighOrderBasisOnMesh::dphi_dxi(dofIndexJ, l, xi)
                  + HighOrderBasisOnMesh::dphi_dxi(dofIndexI, l, xi)*HighOrderBasisOnMesh::dphi_dxi(dofIndexJ, k, xi));
              }
              
              // -----------
              // set kuu
              for (int r = 0; r < 3; r++)
              {  
                for (int s = 0; s < 3; s++)
                {
                  // deps_rs/du_j
                  double depsrs_duj = 1./2*(rightCauchyGreen[s][j]*HighOrderBasisOnMesh::dphi_dxi(dofIndexJ, r, xi)
                    + rightCauchyGreen[r][j]*HighOrderBasisOnMesh::dphi_dxi(dofIndexJ, s, xi));
                 
                  // dp/deps_kl: derivative of \bar{p} w.r.t Green-Lagrange strain tensor E (or eps in paper)
                  double dp_deps_rs = -kappa_*J3*inverseRightCauchyGreen[s][r];
             
                  // d2p/(deps_kl*deps_rs): 2nd derivative of \bar{p} w.r.t Green-Lagrange strain tensor E (or eps in paper), 
                  double factor = 0;
                  for (int c = 0; c < 3; c++)
                  {
                    for (int f = 0; f < 3; f++)
                    {
                      factor += (MathUtility::permutation(k,r,c)*MathUtility::permutation(l,s,f) + MathUtility::permutation(k,s,c)*MathUtility::permutation(l,r,f)) 
                        * rightCauchyGreen[f][c];
                    }
                  }
                  double d2p_deps_klrs = kappa_*J3*inverseRightCauchyGreen[l][k]*inverseRightCauchyGreen[s][r] - kappa_*1./J3*factor;
                  
                  
                  double cuu_klrs = getElasticityEntry(elasticity,k,l,r,s)
                                    + cpp * dp_deps_kl*dp_deps_rs
                                    + cpp * pressureDifference * d2p_deps_klrs;
                                 // + cpp^n * .. terms vanish because d(P(p))/deps = 0
                  double kuu_ij = cuu_klrs*depskl_dui*depsrs_duj + PK2Stress[l][k]*d2epskl_duij;
                  
                  evaluationsKUU[samplingPointIndex](i,j) += kuu_ij;
                  
                }  // s
              }  // r
            
              // -----------
              // set kup 
              // CUP
              double cup_kl = -cpp*dp_deps_kl;
              
              double kup_ij = cup_kl * depskl_dui * dpdpj;
             
              evaluationsKUP[samplingPointIndex](i,j) += kup_ij;
              
            }  // l
          }  // k   
          
          // -----------
          // set kpp 
          double kpp_ij = cpp * dpdpi * dpdpj;
          
          evaluationsKPP[samplingPointIndex](i,j) += kpp_ij;
          
        }  // j
      }  // i
      
      // compute evaluation at current quadrature sampling point for weak form rhs vectors fu, fp
      for (int i = 0; i < nDofsUPerElement; i++)
      {
        int dofIndexI = i;
        // dp_tilde/dp_i where the first p_tilde is the interpolated pressure and p_i is the pressure DOF
        double dpdpi = LowOrderBasisOnMesh::phi(dofIndexI, xi);
        
        // set fu
        // indices of elasticity tensor
        for (int k = 0; k < 3; k++)
        {
          for (int l = 0; l < 3; l++)
          {  
            // deps_kl/du_i
            double depskl_dui = 1./2*(rightCauchyGreen[l][i]*HighOrderBasisOnMesh::dphi_dxi(dofIndexI, k, xi)
                + rightCauchyGreen[k][i]*HighOrderBasisOnMesh::dphi_dxi(dofIndexI, l, xi));
            
            double fu_i = PK2Stress[l][k]*depskl_dui;
            
            evaluationsFU[samplingPointIndex][i] += fu_i;
          }  // l
        }  // k
        
        // set fp
        double fp_i = -cpp * pressureDifference * dpdpi;
        evaluationsFP[samplingPointIndex][i] += fp_i;
      }  // i
      
    }  // function evaluations
  
    // get references to matrices
    Mat &kuu = data_.kuu();
    Mat &kup = data_.kup();
    Mat &kpp = data_.kpp();
    Mat &kppInverse = data_.kppInverse();
    Mat &kupTranspose = data_.kupTranspose();
    Mat &tempKupMatrix = data_.tempKupMatrix();
    Mat &schurComplement = data_.schurComplement();
    Mat &elementStiffness = data_.elementStiffness();
    Mat &stiffnessMatrix = this->data_.stiffnessMatrix();
    Vec &fu = data_.fu();
    Vec &fp = data_.fp();
    Vec &tempKppFp = data_.tempKppFp();
    Vec &tempKupKppFp = data_.tempKupKppFp();
  
    // set entries to 0
    MatZeroEntries(kuu);
    MatZeroEntries(kup);
    MatZeroEntries(kpp);
     
    // perform integration to element stiffness matrices kuu, kup, kpp
    EvaluationsUUType kuuMatrix = QuadratureU::computeIntegral(evaluationsKUU);
    EvaluationsUPType kupMatrix = QuadratureU::computeIntegral(evaluationsKUP);
    EvaluationsPPType kppMatrix = QuadratureP::computeIntegral(evaluationsKPP);
    
    // assign to PETSc types
    kuuMatrix.setPetscMatrix(kuu);
    kupMatrix.setPetscMatrix(kup);
    kppMatrix.setPetscMatrix(kpp);
        
    // perform integration to weak form right hand side
    std::array<double,nDofsUPerElement> fuVector = QuadratureU::computeIntegral(evaluationsFU);
    std::array<double,nDofsPPerElement> fpVector = QuadratureP::computeIntegral(evaluationsFP);
    
    // assign to PETSc types
    PetscUtility::setVector(fuVector,fu);
    PetscUtility::setVector(fpVector,fp);
    
    // -------------
    // by static condensation solve for separate pressure increments
    // explicitly create matrix of schur complement k = kuu - kup * kpp^-1 * kup^T
    MatFactorInfo matFactorInfo;
    IS *permutation;
    ierr = PetscMalloc(sizeof(IS),permutation); CHKERRV(ierr); // allocate permutation struct
    
    // retrive MPI communicator 
    MPI_Comm mpiCommunicator;
    ierr = PetscObjectGetComm((PetscObject)kuu,&mpiCommunicator); CHKERRV(ierr);
    
    // create permutation object
    std::array<int,1> indexSet = {0};
    ierr = ISCreateBlock(mpiCommunicator,nDofsPPerElement,1,indexSet.data(),permutation); CHKERRV(ierr);
    
    // create in-place Cholesky factorization of kpp
    // This changes the state of the matrix to a factored matrix; it cannot be used for example with MatSetValues() unless one first calls MatSetUnfactored(). 
    ierr = MatCholeskyFactor(kpp,*permutation,&matFactorInfo); CHKERRV(ierr);
    
    // compute schur complement kuu - kup*kpp^-1*kup^T
    ierr = MatCreateTranspose(kup, kupTranspose); CHKERRV(ierr);                         // create kup^T
    ierr = MatMatSolve(kpp, kupTranspose, tempKupMatrix); CHKERRV(ierr);                    // solves kpp*t = kup for t = kpp^-1*kup^T
    ierr = MatMatMult(kup, tempKupMatrix, MAT_INITIAL_MATRIX, 1.0, &schurComplement); CHKERRV(ierr);  // s <= kup*t  
    ierr = MatAXPY(schurComplement, -1.0, kuu, SUBSET_NONZERO_PATTERN); CHKERRV(ierr);   // s <= s - kuu
    ierr = MatScale(schurComplement, -1.0); CHKERRV(ierr);                               // s <= -s, i.e. s = -kup*t + kuu
    
    // get the entries of the PETSc matrix in a std::vector
    std::vector<double> schurComplementValues;
    PetscUtility::getMatrixEntries(schurComplement, schurComplementValues);
    
    // assemble local stiffness matrix in global stiffness matrix
    for (int i=0; i<nDofsUPerElement; i++)
    {
      for (int j=0; j<nDofsUPerElement; j++)
      {
        VLOG(2) << "  dof pair (" << i<<","<<j<<"), evaluations: "<<evaluations<<", integrated value: "<<QuadratureU::computeIntegral(evaluations);
        
        double value = schurComplementValues[i*nDofsUPerElement+j];
        ierr = MatSetValue(stiffnessMatrix, dofUNo[i], dofUNo[j], value, ADD_VALUES); CHKERRV(ierr);
      }  // j
    }  // i
    
    // -------------
    // set matching rhs for k
    // F = FU - KUP*KPP^-1*FP
    ierr = MatSolve(kpp, fp, tempKppFp); CHKERRV(ierr);  // solves kpp*t = fp for t = kpp^-1*fp
    ierr = MatMult(kup, tempKppFp, tempKupKppFp); CHKERRV(ierr);  // tempKupKppFp <= kup*kpp^-1*fp  
    ierr = VecAXPY(fu, -1.0, tempKupKppFp); CHKERRV(ierr);   // fu <= fu - tempKupKppFp
    
    std::array<double, nDofsUPerElement> fuValues;
    PetscUtility::getVectorEntries(fu, fuValues);
    
    // assemble local rhs in global rhs
    this->data_.f().setValues(dofUNo, fuValues, ADD_VALUES);
    
  }  // elementNo
  
  // assemble data in rhs vector
  this->data_.f().flushSetValues();
}


// right hand side for compressible mooney rivlin
template<typename BasisOnMeshType, typename MixedQuadratureType>
void FiniteElementMethodStiffnessMatrix<BasisOnMeshType, MixedQuadratureType, Equation::Static::CompressibleMooneyRivlin, Mesh::isDeformable<typename BasisOnMeshType::Mesh>>:: 
void manipulateWeakRhs()
{
  // compute the rhs R-F that fits to the condensed matrix
 
  Vec &rightHandSide = this->data_->rightHandSide();
  Vec &f = this->data_.f().values();
  VecAXPY(rightHandSide, -1.0, f);   // rightHandSide <= rightHandSide - f
}

//! compute the reduced invariants J1 = I1*I3^-1/3, J2 = I2*I3^-2/3, J3=det F
template<typename BasisOnMeshType, typename MixedQuadratureType, typename Term>
std::array<double,3> FiniteElementMethodStiffnessMatrix<BasisOnMeshType, MixedQuadratureType, Term>::
computeReducedInvariants(const std::array<double,3> &invariants, double deformationGradientDeterminant)
{ 
  const double I1 = invariants[0];
  const double I2 = invariants[1];
  const double I3 = invariants[2];  // det C
  
  const double J1 = I1*pow(I3,-1./3);
  const double J2 = I2*pow(I3,-2./3);
  const double J3 = deformationGradientDeterminant; // det F
  
  std::array<double,3> reducedInvariants = {J1, J2, J3};
  
  return reducedInvariants;
}

//! compute the elasticity tensor C = 2*dS/dC. Due to hyperelasticity there are symmetries C_{ijrs} = C_{jirs} and C_{ijrs} = C_{rsij} that leave 21 independent values.
template<typename BasisOnMeshType, typename MixedQuadratureType, typename Term>
std::array<double,21> FiniteElementMethodStiffnessMatrix<BasisOnMeshType, MixedQuadratureType, Term>::
computeElasticityTensor(const std::array<Vec3,3> &rightCauchyGreen, const std::array<Vec3,3> &inverseRightCauchyGreen, const std::array<double,3> invariants, 
                        const std::array<double,3> &reducedInvariants)
{
  // cf. Holzapfel "Nonlinear Solid Mechanics
  
  // the 21 distinct indices (i,j,k,l) of different values of C_{ijkl}
  int indices[21][4] = {
    {0,0,0,0},{0,1,0,0},{0,2,0,0},{1,1,0,0},{1,2,0,0},{2,2,0,0},{0,1,0,1},{0,2,0,1},{1,1,0,1},{1,2,0,1},
    {2,2,0,1},{0,2,0,2},{1,1,0,2},{1,2,0,2},{2,2,0,2},{1,1,1,1},{1,2,1,1},{2,2,1,1},{1,2,1,2},{2,2,1,2},
    {2,2,2,2}
  };
  
  const double I1 = invariants[0];
  const double I2 = invariants[1];
  const double I3 = invariants[2];
  
  const double J1 = reducedInvariants[0];
  const double J2 = reducedInvariants[1];
  const double J3 = reducedInvariants[2];
  
  double rightCauchyGreenDeterminant = I3;
  
  std::array<double,21> elasticity({0});
  // loop over distinct entries in elasticity tensor
  for (int entryNo = 0; entryNo < 21; entryNo++)
  {
    // rename indices of current entry
    const int i = indices[entryNo][0];
    const int j = indices[entryNo][1];
    const int r = indices[entryNo][2];
    const int s = indices[entryNo][3];
    
    // note that rightCauchyGreen and inverseRightCauchyGreen are stored column-major i.e. C_{ij} = rightCauchyGreen[j][i]
    
    // I1''
    double I1dd = 0;
    
    // I2'' = 4*delta_ij*delta_rs - 2*(delta_ir*delta_js + delta_is*delta_jr)
    double I2dd = 4*(i==j)*(r==s) - 2*((i==r)*(j==s) + (i==s)*(j==r));
    
    // I3''
    double I3dd = 0;
    for (int c = 0; c < 3; c++)
    {
      for (int f = 0; f < 3; f++)
      {
        I3dd += (MathUtility::permutation(i,r,c)*MathUtility::permutation(j,s,f)
                   + MathUtility::permutation(i,s,c)*MathUtility::permutation(j,r,f)
                   + MathUtility::permutation(j,r,c)*MathUtility::permutation(i,s,f)
                   + MathUtility::permutation(j,s,c)*MathUtility::permutation(i,r,f)) * rightCauchyGreen[f][c];
      }
    }
    
    // I1' = 2*delta_ij
    double I1dij = 2*(i==j);
    double I1drs = 2*(r==s);
    
    // I2' 
    double I2dij = 2*I1*(i==j) - (rightCauchyGreen[j][i] + rightCauchyGreen[i][j]);
    double I2drs = 2*I1*(r==s) - (rightCauchyGreen[s][r] + rightCauchyGreen[r][s]);
    
    // I3'
    double I3dij = (inverseRightCauchyGreen[j][i] + inverseRightCauchyGreen[i][j]) * rightCauchyGreenDeterminant;
    double I3drs = (inverseRightCauchyGreen[s][r] + inverseRightCauchyGreen[r][s]) * rightCauchyGreenDeterminant;
    
    // J1''
    double J1dd = pow(I3,-1./3)*I1dd - 1./3*pow(I3,-4./3) * (I1dij*I3drs + I3dij*I1drs + I1*I3dd)
     + 4./9.*(I1*pow(I3,-7./3))*I3dij*I3drs;
     
    // J2''
    double J2dd = pow(I3,-2./3)*I2dd - 2./3*pow(I3,-5./3) * (I2dij*I3drs + I3dij*I2drs + I2*I3dd)
     + 10./9.*(I2*pow(I3,-8./3))*I3dij*I3drs;
    
    // J3''
    double J3dd = -1./4*pow(I3,-3./2)*I3dij*I3drs + 1./2*pow(I3,-1./2)*I3dd;
     
    // J1'
    double J1d = pow(I3,-1./3)*I1dij - 1./3*(I1*pow(I3,-4./3))*I3dij;
    
    // J2'
    double J2d = pow(I3,-2./3)*I2dij - 2./3*(I2*pow(I3,-5./3))*I3dij;
    
    // J3'
    double J3dij = 1./2*pow(I3,-1./2)*I3dij;
    double J3drs = 1./2*pow(I3,-1./2)*I3drs;
    
    elasticity[entryNo] = c1_*J1dd + c2_*J2dd + kappa_*(J3dij*J3drs + (J3 - 1)*J3dd);
    
  }
  
  return elasticity;
}

//! compute 2nd Piola-Kirchhoff stress tensor S = 2*sym(dPsi/dC)template<typename BasisOnMeshType, typename MixedQuadratureType, typename Term>
std::array<Vec3,3> FiniteElementMethodStiffnessMatrix<BasisOnMeshType, MixedQuadratureType, Term>::
computePK2Stress(const std::array<Vec3,3> &rightCauchyGreen, const std::array<Vec3,3> &inverseRightCauchyGreen, const std::array<double,3> invariants, 
                 const std::array<double,3> &reducedInvariants)
{
  const double I1 = invariants[0];
  const double I2 = invariants[1];
  const double I3 = invariants[2];
  
  const double J1 = reducedInvariants[0];
  const double J2 = reducedInvariants[1];
  const double J3 = reducedInvariants[2];
  
  double rightCauchyGreenDeterminant = I3;
   
  std::array<Vec3,3> pK2Stress({0});
 
  // row index
  for (int i=0; i<3; i++)
  {
    // column index
    for (int j=0; j<3; j++)
    { 
      // I1' = 2*delta_ij
      double I1dij = 2*(i==j);
      
      // I2' 
      double I2dij = 2*I1*(i==j) - (rightCauchyGreen[j][i] + rightCauchyGreen[i][j]);
      
      // I3'
      double I3dij = (inverseRightCauchyGreen[j][i] + inverseRightCauchyGreen[i][j]) * rightCauchyGreenDeterminant;
     
      // J1'
      double J1dij = pow(I3,-1./3)*I1dij - 1./3*(I1*pow(I3,-4./3))*I3dij;
      
      // J2'
      double J2dij = pow(I3,-2./3)*I2dij - 2./3*(I2*pow(I3,-5./3))*I3dij;

      // J3'
      double J3dij = 1./2*pow(I3,-1./2)*I3dij;
    
      pK2Stress[j][i] = c1_*J1dij + c2_*J2dij + kappa_*(J3 - 1)*J3dij;
    }
  }
}

//! return the entry klrs of the elasticity tensor
template<typename BasisOnMeshType, typename MixedQuadratureType, typename Term>
int FiniteElementMethodStiffnessMatrix<BasisOnMeshType, MixedQuadratureType, Term>::
getElasticityEntryNo(int k, int l, int r, int s)
{
  // this method was tested outside of this codebase and is correct
  // the inverse mapping is given by the following array:
  //int indices[21][4] = {
  //  {0,0,0,0},{0,1,0,0},{0,2,0,0},{1,1,0,0},{1,2,0,0},{2,2,0,0},{0,1,0,1},{0,2,0,1},{1,1,0,1},{1,2,0,1},
  //  {2,2,0,1},{0,2,0,2},{1,1,0,2},{1,2,0,2},{2,2,0,2},{1,1,1,1},{1,2,1,1},{2,2,1,1},{1,2,1,2},{2,2,1,2},
  //  {2,2,2,2}
  //};
  switch(k)
  {
  case 0: // k
    switch(l)
    {
    case 0: // l
      switch(r)
      {
      case 0: // 000
        switch(s)
        {
        case 0:
          return 0;
          break;
        case 1:
          return 1;
          break;
        case 2:
          return 2;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 1: // 001
        switch(s)
        {
        case 0:
          return 1;
          break;
        case 1:
          return 3;
          break;
        case 2:
          return 4;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 2: // 002
        switch(s)
        {
        case 0:
          return 2;
          break;
        case 1:
          return 4;
          break;
        case 2:
          return 5;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      default:
        LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
        break;
      }  // r
      break;
    case 1: // l
      switch(r)
      {
      case 0: // 010
        switch(s)
        {
        case 0:
          return 1;
          break;
        case 1:
          return 6;
          break;
        case 2:
          return 7;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 1: // 011
        switch(s)
        {
        case 0:
          return 6;
          break;
        case 1:
          return 8;
          break;
        case 2:
          return 9;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 2: // 012
        switch(s)
        {
        case 0:
          return 7;
          break;
        case 1:
          return 9;
          break;
        case 2:
          return 10;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      default:
        LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
        break;
      }  // r
      break;
    case 2: // l
      switch(r)
      {
      case 0: // 020
        switch(s)
        {
        case 0:
          return 2;
          break;
        case 1:
          return 7;
          break;
        case 2:
          return 11;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 1: // 021
        switch(s)
        {
        case 0:
          return 7;
          break;
        case 1:
          return 12;
          break;
        case 2:
          return 13;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 2: // 022
        switch(s)
        {
        case 0:
          return 11;
          break;
        case 1:
          return 13;
          break;
        case 2:
          return 14;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      default:
        LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
        break;
      }  // r
      break;
    default:
      LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
      break;
    }  // l
    break;
  case 1: // k
    switch(l)
    {
    case 0: // l
      switch(r)
      {
      case 0: // 100
        switch(s)
        {
        case 0:
          return 1;
          break;
        case 1:
          return 6;
          break;
        case 2:
          return 7;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 1: // 101
        switch(s)
        {
        case 0:
          return 6;
          break;
        case 1:
          return 8;
          break;
        case 2:
          return 9;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 2: // 102
        switch(s)
        {
        case 0:
          return 7;
          break;
        case 1:
          return 9;
          break;
        case 2:
          return 10;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      default:
        LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
        break;
      }  // r
      break;
    case 1: // l
      switch(r)
      {
      case 0: // 110
        switch(s)
        {
        case 0:
          return 3;
          break;
        case 1:
          return 8;
          break;
        case 2:
          return 12;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 1: // 111
        switch(s)
        {
        case 0:
          return 8;
          break;
        case 1:
          return 15;
          break;
        case 2:
          return 16;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 2: // 112
        switch(s)
        {
        case 0:
          return 12;
          break;
        case 1:
          return 16;
          break;
        case 2:
          return 17;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      default:
        LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
        break;
      }  // r
      break;
    case 2: // l
      switch(r)
      {
      case 0: // 120
        switch(s)
        {
        case 0:
          return 4;
          break;
        case 1:
          return 9;
          break;
        case 2:
          return 13;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 1: // 121
        switch(s)
        {
        case 0:
          return 9;
          break;
        case 1:
          return 16;
          break;
        case 2:
          return 18;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 2: // 122
        switch(s)
        {
        case 0:
          return 13;
          break;
        case 1:
          return 18;
          break;
        case 2:
          return 19;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      default:
        LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
        break;
      }  // r
      break;
    default:
      LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
      break;
    }  // l
    break;
  case 2: // k
    switch(l)
    {
    case 0: // l
      switch(r)
      {
      case 0: // 200
        switch(s)
        {
        case 0:
          return 2;
          break;
        case 1:
          return 7;
          break;
        case 2:
          return 11;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 1: // 201
        switch(s)
        {
        case 0:
          return 7;
          break;
        case 1:
          return 12;
          break;
        case 2:
          return 13;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 2: // 202
        switch(s)
        {
        case 0:
          return 11;
          break;
        case 1:
          return 13;
          break;
        case 2:
          return 14;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      default:
        LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
        break;
      }  // r
      break;
    case 1: // l
      switch(r)
      {
      case 0: // 210
        switch(s)
        {
        case 0:
          return 4;
          break;
        case 1:
          return 9;
          break;
        case 2:
          return 13;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 1: // 211
        switch(s)
        {
        case 0:
          return 9;
          break;
        case 1:
          return 16;
          break;
        case 2:
          return 18;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 2: // 212
        switch(s)
        {
        case 0:
          return 13;
          break;
        case 1:
          return 18;
          break;
        case 2:
          return 19;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      default:
        LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
        break;
      }  // r
      break;
    case 2: // l
      switch(r)
      {
      case 0: // 220
        switch(s)
        {
        case 0:
          return 5;
          break;
        case 1:
          return 10;
          break;
        case 2:
          return 14;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 1: // 221
        switch(s)
        {
        case 0:
          return 10;
          break;
        case 1:
          return 17;
          break;
        case 2:
          return 19;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      case 2: // 222
        switch(s)
        {
        case 0:
          return 14;
          break;
        case 1:
          return 19;
          break;
        case 2:
          return 20;
          break;
        default:
          LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
          break;
        }  // s
        break;
      default:
        LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
        break;
      }  // r
      break;
    default:
      LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
      break;
    }  // l
    break;
  default:
    LOG(ERROR) << "Wrong index (" << k << "," << l << "," << r << "," << s << ") to elasticity tensor!";
    break;
  }  // k
  return 0;
}


//! return the entry klrs of the elasticity tensor
template<typename BasisOnMeshType, typename MixedQuadratureType, typename Term>
double FiniteElementMethodStiffnessMatrix<BasisOnMeshType, MixedQuadratureType, Term>::
getElasticityEntry(std::array<double, 21> &elasticity, int k, int l, int r, int s)
{
  return elasticity[getElasticityEntryNo(k,l,r,s)];
}
*/
};    // namespace