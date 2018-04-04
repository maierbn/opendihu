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
  const int nUnknowsPerElement = nDofsPerElement*3;    // 3 directions for displacements per dof
  
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
      Tensor2 jacobianMaterial = BasisOnMeshType::computeJacobian(geometryReferenceValues, xi);
      Tensor2 inverseJacobianMaterial = MathUtility::computeInverse(jacobianMaterial);
      // jacobianMaterial[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx
      
      // F
      Tensor2 deformationGradient = this->computeDeformationGradient(displacementValues, jacobianMaterial, xi);
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
      Tensor4 elasticity = this->computeElasticityTensor(artificialPressure, artificialPressureTilde, rightCauchyGreen, inverseRightCauchyGreen, fictitiousPK2Stress, pk2StressIsochoric, deformationGradientDeterminant);
      
      // get evaluations of integrand at xi for all (i,j)-dof pairs, integrand is defined in another class
      //evaluationsArray[samplingPointIndex] = 
         
      std::array<Vec3,nDofsPerElement> gradPhi = mesh->getGradPhi(xi);
      // (column-major storage) gradPhi[M][a] = dphi_M / dxi_a
      
      // loop over pairs of basis functions and evaluation integrand at xi
      for (int aDof = 0; aDof < nDofsPerElement; aDof++)           // index over dofs, each dof has D components
      {
        for (int aComponent = 0; aComponent < D; aComponent++)     // lower-case a in derivation
        {
          for (int bInternal = 0; bInternal < D; bInternal++)     // capital B in derivation
          {
            // compute index of degree of freedom and component (matrix row index)
            const int i = nDofsPerElement*D*aDof + D*aComponent + bInternal;
           
            for (int bDof = 0; bDof < nDofsPerElement; bDof++)
            {
              for (int bComponent = 0; bComponent < D; bComponent++)     // lower-case b in derivation
              {
                for (int dInternal = 0; dInternal < D; dInternal++)     // capital D in derivation
                {
                  // compute index of degree of freedom and component (index of entry in vector of unknows / matrix column index)
                  const int j = nDofsPerElement*D*bDof + D*bComponent + dInternal;
              
                  // k_ij = int dphi_a/dX_B * dphi_b/dX_D * k_abBD 
                  
                  // gradPhi[columnIdx][rowIdx] = dphi_(dof columnIdx)/dxi_rowIdx (independent of displacement component)
                  // jacobianMaterial[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx (independent of displacement component)
                  
                  
                  double integrand = gradPhi[aDof][bInternal] * gradPhi[bDof][dInternal];
                  //evaluations[i][j] = integrand;
                }
              }
            }
          }
        }
      }
      
    }  // function evaluations
    
    // perform integration and add to entry of stiffness matrix
    for (int a=0; a<nDofsUPerElement; a++)
    {
      for (int b=0; b<nDofsUPerElement; b++)
      {
        // extract evaluations for current (a,b) dof-pair
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