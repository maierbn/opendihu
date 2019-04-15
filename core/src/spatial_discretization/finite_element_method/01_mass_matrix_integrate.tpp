#include "spatial_discretization/finite_element_method/01_matrix.h"

#include <Python.h>
#include <memory>
#include <petscksp.h>
#include <petscsys.h>

#include "quadrature/tensor_product.h"
#include "spatial_discretization/finite_element_method/integrand/integrand_mass_matrix.h"

namespace SpatialDiscretization
{

// 1D,2D,3D mass matrix of Deformable mesh
template<typename FunctionSpaceType,typename QuadratureType,typename Term,typename Dummy0,typename Dummy1,typename Dummy2>
void FiniteElementMethodMatrix<FunctionSpaceType,QuadratureType,Term,Dummy0,Dummy1,Dummy2>::
setMassMatrix()
{
  // check if matrix discretization matrix exists
  if (!this->data_.massMatrix())
  {
    this->data_.initializeMassMatrix();
  }

  const int D = FunctionSpaceType::dim();
  LOG(TRACE) << "setMassMatrix " << D << "D using integration, FunctionSpaceType: " << typeid(FunctionSpaceType).name() << ", QuadratureType: " << typeid(QuadratureType).name();

  // massMatrix * f_strong = rhs_weak
  // row of massMatrix: contributions to a single entry in rhs_weak

  // define shortcuts for integrator and basis
  typedef Quadrature::TensorProduct<D,QuadratureType> QuadratureDD;
  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();
  typedef MathUtility::Matrix<nDofsPerElement,nDofsPerElement> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureDD::numberEvaluations()
          > EvaluationsArrayType;     // evaluations[nGP^D][nDofs][nDofs]

  // initialize variables
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> massMatrix = this->data_.massMatrix();

  std::shared_ptr<FunctionSpaceType> functionSpace = std::static_pointer_cast<FunctionSpaceType>(this->data_.functionSpace());
  functionSpace->geometryField().setRepresentationGlobal();
  functionSpace->geometryField().startGhostManipulation();   // ensure that local ghost values of geometry field are set

  // initialize values to zero
  // loop over elements
  for (element_no_t elementNo = 0; elementNo < functionSpace->nElementsLocal(); elementNo++)
  {
    std::array<dof_no_t,nDofsPerElement> dofNosLocal = functionSpace->getElementDofNosLocal(elementNo);

    for (int i=0; i<nDofsPerElement; i++)
    {
      for (int j=0; j<nDofsPerElement; j++)
      {
        massMatrix->setValue(dofNosLocal[i], dofNosLocal[j], 0, INSERT_VALUES);
      }
    }
  }
  massMatrix->assembly(MAT_FLUSH_ASSEMBLY);

  // setup arrays used for integration
  std::array<std::array<double,D>, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();
  EvaluationsArrayType evaluationsArray{};

  LOG(DEBUG) << "1D integration with " << QuadratureType::numberEvaluations() << " evaluations";
  LOG(DEBUG) << D << "D integration with " << QuadratureDD::numberEvaluations() << " evaluations";

  // set entries in massMatrix
  // loop over elements
  for (element_no_t elementNo = 0; elementNo < functionSpace->nElementsLocal(); elementNo++)
  {
    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNosLocal = functionSpace->getElementDofNosLocal(elementNo);

    // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
    std::array<Vec3,FunctionSpaceType::nDofsPerElement()> geometry;
    functionSpace->getElementGeometry(elementNo, geometry);

    // compute integral
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // evaluate function to integrate at samplingPoints[i*2], write value to evaluations[i]
      std::array<double,D> xi = samplingPoints[samplingPointIndex];

      // compute the 3xD jacobian of the parameter space to world space mapping
      auto jacobian = FunctionSpaceType::computeJacobian(geometry, xi);

      // get evaluations of integrand which is defined in another class
      evaluationsArray[samplingPointIndex] = IntegrandMassMatrix<D,EvaluationsType,FunctionSpaceType,Term>::evaluateIntegrand(jacobian,xi);

    }  // function evaluations

    // integrate all values for the (i,j) dof pairs at once
    EvaluationsType integratedValues = QuadratureDD::computeIntegral(evaluationsArray);

    // perform integration and add to entry in rhs vector
    for (int i=0; i<nDofsPerElement; i++)
    {
      for (int j=0; j<nDofsPerElement; j++)
      {
        // integrate value and set entry in discretization matrix
        double integratedValue = integratedValues(i,j);

        massMatrix->setValue(dofNosLocal[i], dofNosLocal[j], integratedValue, ADD_VALUES);
      }  // j
    }  // i
  }  // elementNo

  // merge local changes in parallel and assemble the matrix (MatAssemblyBegin, MatAssemblyEnd)
  massMatrix->assembly(MAT_FINAL_ASSEMBLY);
}

}  // namespace
