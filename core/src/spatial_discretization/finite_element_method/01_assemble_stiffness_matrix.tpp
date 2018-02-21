#include "spatial_discretization/finite_element_method/01_assemble_stiffness_matrix.h"

#include <Python.h>
#include <memory>
#include <vector>
#include <petscsys.h>

#include "quadrature/tensor_product.h"
#include "basis_on_mesh/05_basis_on_mesh.h"
#include "spatial_discretization/finite_element_method/01_integrand_stiffness_matrix.h"

namespace SpatialDiscretization
{
  
// 1D,2D,3D stiffness matrix of Deformable mesh
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void AssembleStiffnessMatrix<BasisOnMeshType, QuadratureType, Term>::
setStiffnessMatrix()
{
  const int D = BasisOnMeshType::dim();
  LOG(TRACE)<<"setStiffnessMatrix " << D << "D";
 
  // get prefactor value
  const double prefactor = PythonUtility::getOptionDouble(this->specificSettings_, "prefactor", 1.0);
  
  // define shortcuts for integrator and basis
  typedef Quadrature::TensorProduct<D,QuadratureType> QuadratureDD;
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  typedef std::array<std::array<double, nDofsPerElement>, nDofsPerElement> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureDD::numberEvaluations()
          > EvaluationsArrayType;     // evaluations[nGP^D][nDofs][nDofs]
  
  // initialize variables
  PetscErrorCode ierr;
  Mat &stiffnessMatrix = this->data_.stiffnessMatrix();
  
  std::shared_ptr<BasisOnMeshType> mesh = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh());
  
  // initialize values to zero
  int cntr = 1;
  // loop over elements 
  for (element_no_t elementNo = 0; elementNo < mesh->nElements(); elementNo++)
  {
    auto dof = mesh->getElementDofNos(elementNo);
    
    for (int i=0; i<nDofsPerElement; i++)
    {
      for (int j=0; j<nDofsPerElement; j++)
      {
        VLOG(3) << " initialize stiffnessMatrix entry ( " << dof[i] << "," << dof[j] << ") (no. " << cntr++ << ")";
        ierr = MatSetValue(stiffnessMatrix, dof[i], dof[j], 0, INSERT_VALUES); CHKERRV(ierr);
      }
    }
  }
  
  // setup arrays used for integration
  std::array<std::array<double,D>, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();
  EvaluationsArrayType evaluationsArray;
  
  LOG(DEBUG) << "1D integration with " << QuadratureType::numberEvaluations() << " evaluations";
  LOG(DEBUG) << D << "D integration with " << QuadratureDD::numberEvaluations() << " evaluations";
#ifdef DEBUG
  LOG(DEBUG) << "SAMPLING POINTS: ";
  for  (auto value : samplingPoints)
    LOG(DEBUG) << "   " << value;
#endif
  
  // fill entries in stiffness matrix
  // loop over elements 
  for (int elementNo = 0; elementNo < mesh->nElements(); elementNo++)
  {
    // get indices of element-local dofs
    auto dof = mesh->getElementDofNos(elementNo);
    
    VLOG(2) << "element " << elementNo;
    
    // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
    std::array<Vec3,BasisOnMeshType::nDofsPerElement()> geometry;
    mesh->getElementGeometry(elementNo, geometry);
    
    // compute integral
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // evaluate function to integrate at samplingPoint
      std::array<double,D> xi = samplingPoints[samplingPointIndex];
      
      // compute the 3xD jacobian of the parameter space to world space mapping
      auto jacobian = BasisOnMeshType::computeJacobian(geometry, xi);
      
      VLOG(2) << "samplingPointIndex="<<samplingPointIndex<<", xi="<<xi<<", geometry: "<<geometry<<", jac: " <<jacobian;
      
      // get evaluations of integrand at xi for all (i,j)-dof pairs, integrand is defined in another class
      evaluationsArray[samplingPointIndex] 
        = IntegrandStiffnessMatrix<D,EvaluationsType,BasisOnMeshType,Term>::evaluateIntegrand(mesh,jacobian,xi);
      
      
    }  // function evaluations
    
    // perform integration and add to entry of stiffness matrix
    for (int i=0; i<nDofsPerElement; i++)
    {
      for (int j=0; j<nDofsPerElement; j++)
      {
        // extract evaluations for current (i,j) dof-pair
        std::array<double,QuadratureDD::numberEvaluations()> evaluations;
        for (int k=0; k<QuadratureDD::numberEvaluations(); k++)
          evaluations[k] = evaluationsArray[k][i][j];
        
        VLOG(2) << "  dof pair (" << i<<","<<j<<"), evaluations: "<<evaluations<<", integrated value: "<<QuadratureDD::computeIntegral(evaluations);
        
        // integrate value and set entry in stiffness matrix
        double value = -prefactor * QuadratureDD::computeIntegral(evaluations);
        ierr = MatSetValue(stiffnessMatrix, dof[i], dof[j], value, ADD_VALUES); CHKERRV(ierr);
      }  // j
    }  // i
  }  // elementNo
}

};
