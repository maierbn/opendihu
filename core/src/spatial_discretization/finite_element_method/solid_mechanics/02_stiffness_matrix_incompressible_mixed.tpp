#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_incompressible.h"

#include <Python.h>  // has to be the first included header
#include <memory>
#include <petscksp.h>
#include <petscsys.h>

#include "semt/Semt.h"
#include "semt/Shortcuts.h"

namespace SpatialDiscretization
{
  
// set stiffness matrix for a u-p mixed formulation in which the pressure is not condensed out
template<typename LowOrderBasisOnMeshType, typename HighOrderBasisOnMeshType, typename MixedQuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,
  MixedQuadratureType, 
  Term,
  std::enable_if_t<LowOrderBasisOnMeshType::BasisFunction::isNodalBased, typename HighOrderBasisOnMeshType::Mesh>,
  Equation::isIncompressible<Term>
>::
setStiffnessMatrix()
{
  LOG(TRACE)<<"setStiffnessMatrix for solid mechanics";
/*
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
  
  typedef std::array<std::array<double, nDofsUPerElement>, nDofsUPerElement> EvaluationsUType;
  typedef std::array<
            EvaluationsUType,
            QuadratureU::numberEvaluations()
          > EvaluationsUArrayType;     // evaluations[nGP^D][nDofs][nDofs]
  typedef std::array<std::array<double, nDofsPPerElement>, nDofsPPerElement> EvaluationsPType;
  typedef std::array<
            EvaluationsPType,
            QuadratureP::numberEvaluations()
          > EvaluationsPArrayType;     // evaluations[nGP^D][nDofs][nDofs]
  
  // setup arrays used for integration
  std::array<std::array<double,D>, QuadratureU::numberEvaluations()> samplingPointsU = QuadratureU::samplingPoints();
  std::array<std::array<double,D>, QuadratureP::numberEvaluations()> samplingPointsP = QuadratureP::samplingPoints();
  EvaluationsUArrayType evaluationsUArray;
  //EvaluationsPArrayType evaluationsPArray;
  
  typedef std::array<Vec3,BasisOnMeshType::dim()> Tensor2;
  typedef std::array<double,21> Tensor4;   // data type for 4th order elasticity tensor with due to symmetry has 21 independent components
  
  // loop over elements 
  for (int elementNo = 0; elementNo < nElements; elementNo++)
  {
    // get indices of element-local dofs
    auto dofUNo = basisOnMeshU->getElementDofNos(elementNo);
    auto dofPNo = basisOnMeshP->getElementDofNos(elementNo);
    
    // get geometry field of current configuration of meshU
    std::array<Vec3,HighOrderBasisOnMesh::nDofsPerElement()> geometryCurrent;
    basisOnMeshU->getElementGeometry(elementNo, geometryCurrent);
        
    // get geometry field of reference configuration
    std::array<Vec3,HighOrderBasisOnMesh::nDofsPerElement()> geometryReference;
    this->data_.geometryReference().getElementValues(elementNo, geometryReference);
    
    // get displacement field
    std::array<Vec3,HighOrderBasisOnMesh::nDofsPerElement()> displacement;
    this->data_.displacement().getElementValues(elementNo, displacement);
    
    // get separately interpolated pressure field
    std::array<double,LowOrderBasisOnMesh::nDofsPerElement()> separatePressure;
    this->data_.pressure().getElementValues(elementNo, separatePressure);
    
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
      
      // Pk2 stress tensor S = 
      Tensor2 PK2Stress = this->computePK2Stress(rightCauchyGreen, inverseRightCauchyGreen, invariants);
      
      // elasticity tensor C_{ijkl}
      Tensor4 elasticity = this->computeElasticityTensor(rightCauchyGreen, inverseRightCauchyGreen, invariants);
      
      // p = -1/(3J)C:S
      double pressureFromDisplacements = this->computePressureFromDisplacements(deformationGradientDeterminant, rightCauchyGreen, PK2Stress);
      pressureFromDisplacements++;  // avoid warning unused variable
      VLOG(2) << "pressureFromDisplacements: " << pressureFromDisplacements;
      VLOG(2) << "samplingPointIndex="<<samplingPointIndex<<", xi="<<xi;
      VLOG(2) << "geometryCurrent="<<geometryCurrent<<", geometryReference="<<geometryReference;
      
      // get evaluations of integrand at xi for all (i,j)-dof pairs, integrand is defined in another class
      //evaluationsArray[samplingPointIndex] = 
      
    }  // function evaluations
    
    // perform integration and add to entry of stiffness matrix
    for (int i=0; i<nDofsUPerElement; i++)
    {
      for (int j=0; j<nDofsUPerElement; j++)
      {
        // extract evaluations for current (i,j) dof-pair
        std::array<double,QuadratureU::numberEvaluations()> evaluations;
        for (int k=0; k<QuadratureU::numberEvaluations(); k++)
          evaluations[k] = evaluationsUArray[k][i][j];
        
        VLOG(2) << "  dof pair (" << i<<","<<j<<"), evaluations: "<<evaluations<<", integrated value: "<<QuadratureU::computeIntegral(evaluations);
        
        // integrate value and set entry in stiffness matrix
        //double value = QuadratureU::computeIntegral(evaluations);
      }  // j
    }  // i
  }  // elementNo
  */
}

};    // namespace