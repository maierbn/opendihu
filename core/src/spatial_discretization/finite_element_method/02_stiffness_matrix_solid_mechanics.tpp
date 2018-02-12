#include "spatial_discretization/finite_element_method/02_stiffness_matrix.h"

#include <Python.h>  // has to be the first included header
#include <memory>
#include <petscksp.h>
#include <petscsys.h>

namespace SpatialDiscretization
{
  
template<typename BasisOnMeshType, typename MixedIntegratorType, typename Term>
void FiniteElementMethodStiffnessMatrix<BasisOnMeshType, MixedIntegratorType, Term, Mesh::isDeformable<typename BasisOnMeshType::Mesh>, Equation::isSolidMechanics<Term>>:: 
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
  typedef Integrator::TensorProduct<D,typename MixedIntegratorType::HighOrderIntegrator> IntegratorU;
  typedef Integrator::TensorProduct<D,typename MixedIntegratorType::LowOrderIntegrator> IntegratorP;
  
  typedef std::array<std::array<double, nDofsUPerElement>, nDofsUPerElement> EvaluationsUType;
  typedef std::array<
            EvaluationsUType,
            IntegratorU::numberEvaluations()
          > EvaluationsUArrayType;     // evaluations[nGP^D][nDofs][nDofs]
  typedef std::array<std::array<double, nDofsPPerElement>, nDofsPPerElement> EvaluationsPType;
  typedef std::array<
            EvaluationsPType,
            IntegratorP::numberEvaluations()
          > EvaluationsPArrayType;     // evaluations[nGP^D][nDofs][nDofs]
  
  // setup arrays used for integration
  std::array<std::array<double,D>, IntegratorU::numberEvaluations()> samplingPointsU = IntegratorU::samplingPoints();
  std::array<std::array<double,D>, IntegratorP::numberEvaluations()> samplingPointsP = IntegratorP::samplingPoints();
  EvaluationsUArrayType evaluationsUArray;
  EvaluationsPArrayType evaluationsPArray;
  
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
    this->data_.geometryReference().template getElementValues<D>(elementNo, geometryReference);
    
    // get displacement field
    std::array<Vec3,HighOrderBasisOnMesh::nDofsPerElement()> displacement;
    this->data_.displacement().template getElementValues<D>(elementNo, displacement);
    
    //TODO: generate geometryReference in initializations
    // loop over integration points (e.g. gauss points) for displacement field
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPointsU.size(); samplingPointIndex++)
    {
      // get parameter values of current sampling point
      std::array<double,D> xi = samplingPointsU[samplingPointIndex];
      
      // compute the 3xD jacobian of the parameter space to world space mapping
      auto jacobian = HighOrderBasisOnMesh::computeJacobian(geometryCurrent, xi);
      
      auto deformationGradient = computeDeformationGradient(geometryReference, displacement, jacobian, xi);
      
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
        std::array<double,IntegratorU::numberEvaluations()> evaluations;
        for (int k=0; k<IntegratorU::numberEvaluations(); k++)
          evaluations[k] = evaluationsUArray[k][i][j];
        
        VLOG(2) << "  dof pair (" << i<<","<<j<<"), evaluations: "<<evaluations<<", integrated value: "<<IntegratorU::integrate(evaluations);
        
        // integrate value and set entry in stiffness matrix
        //double value = IntegratorU::integrate(evaluations);
      }  // j
    }  // i
  }  // elementNo
}

template<typename BasisOnMeshType, typename MixedIntegratorType, typename Term>
std::array<Vec3,BasisOnMeshType::dim()> FiniteElementMethodStiffnessMatrix<BasisOnMeshType, MixedIntegratorType, Term, Mesh::isDeformable<typename BasisOnMeshType::Mesh>, Equation::isSolidMechanics<Term>>:: 
computeDeformationGradient(std::array<Vec3,BasisOnMeshType::HighOrderBasisOnMesh::nDofsPerElement()> &geometryReference,
                           std::array<Vec3,BasisOnMeshType::HighOrderBasisOnMesh::nDofsPerElement()> &displacement, 
                           std::array<Vec3,BasisOnMeshType::dim()> &jacobian, 
                           std::array<double, BasisOnMeshType::dim()> xi)
{
}


};    // namespace