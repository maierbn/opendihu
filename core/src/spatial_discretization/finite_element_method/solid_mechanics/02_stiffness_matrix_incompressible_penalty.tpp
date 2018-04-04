#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_incompressible.h"

#include <Python.h>  // has to be the first included header
#include <memory>
#include <petscksp.h>
#include <petscsys.h>

#include "semt/Semt.h"
#include "semt/Shortcuts.h"
#include "control/types.h"

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
  
  const int D = BasisOnMeshType::dim();
  const int nDofs = mesh->nDofs();
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  const int nElements = mesh->nElements();
  
  // define shortcuts for integrator and basis
  typedef Quadrature::TensorProduct<D,QuadratureType> QuadratureDD;
  
  typedef std::array<std::array<double, nDofsPerElement>, nDofsPerElement> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureType::numberEvaluations()
          > EvaluationsArrayType;     // evaluations[nGP^D][nDofs][nDofs]
  
  // setup arrays used for integration
  std::array<std::array<double,D>, Quadrature::numberEvaluations()> samplingPoints = QuadratureType::samplingPoints();
  EvaluationsArrayType evaluationsArray;
  
  typedef std::array<Vec3,BasisOnMeshType::dim()> Tensor2;
  typedef std::array<double,21> Tensor4;   // data type for 4th order elasticity tensor with due to symmetry has 21 independent components
  
  
  // loop over elements 
  for (int elementNo = 0; elementNo < nElements; elementNo++)
  {
    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNo = mesh->getElementDofNos(elementNo);
    
    // get geometry field of current configuration of mesh
    std::array<Vec3,nDofsPerElement> geometryCurrentValues;
    mesh->getElementGeometry(elementNo, geometryCurrentValues);
        
    // get geometry field of reference configuration
    std::array<Vec3,nDofsPerElement> geometryReferenceValues;
    this->data_.geometryReference().getElementValues(elementNo, geometryReferenceValues);
    
    // get displacement field
    std::array<Vec3,nDofsPerElement> displacementValues;
    this->data_.displacement().getElementValues(elementNo, displacementValues);
    
    // loop over integration points (e.g. gauss points) for displacement field
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // get parameter values of current sampling point
      std::array<double,D> xi = samplingPoints[samplingPointIndex];
      
      // compute the 3xD jacobian of the parameter space to world space mapping
      Tensor2 jacobian = BasisOnMeshType::computeJacobian(geometryCurrentValues, xi);
      
      // F
      Tensor2 deformationGradient = this->computeDeformationGradient(displacementValues, jacobian, xi);
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
      Tensor2 PK2Stress = this->computePK2Stress(artificialPressure, rightCauchyGreen, inverseRightCauchyGreen, reducedInvariants, deformationGradientDeterminant, fictitiousPK2Stress);
      
      // elasticity tensor C_{ijkl}
      Tensor4 elasticity = this->computeElasticityTensor(artificialPressure, artificialPressureTilde, rightCauchyGreen, inverseRightCauchyGreen, fictitiousPK2Stress, deformationGradientDeterminant);
      
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
  
}

};    // namespace