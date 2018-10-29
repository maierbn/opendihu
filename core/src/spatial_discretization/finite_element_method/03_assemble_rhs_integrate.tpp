#include "spatial_discretization/finite_element_method/03_assemble_rhs.h"

#include <Python.h>
#include <memory>
#include <vector>
#include <petscsys.h>

#include "quadrature/tensor_product.h"
#include "function_space/function_space.h"
#include "spatial_discretization/finite_element_method/integrand/integrand_mass_matrix.h"
#include "field_variable/field_variable.h"

namespace SpatialDiscretization
{

// 1D,2D,3D rhs vector of Deformable mesh
template<typename FunctionSpaceType, typename QuadratureType, typename Term, typename Dummy>
void AssembleRightHandSide<FunctionSpaceType, QuadratureType, Term, Dummy>::
multiplyRightHandSideWithMassMatrix()
{
  const int D = FunctionSpaceType::dim();
  LOG(TRACE) << "multiplyRightHandSideWithMassMatrix " << D << "D";

  // define shortcuts for integrator and basis
  typedef Quadrature::TensorProduct<D,QuadratureType> QuadratureDD;
  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();
  typedef MathUtility::Matrix<nDofsPerElement,nDofsPerElement> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureDD::numberEvaluations()
          > EvaluationsArrayType;    // evaluations[nGP^D][nDofs][nDofs]

  // initialize variables
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> rightHandSide = this->data_.rightHandSide();

  std::shared_ptr<FunctionSpaceType> functionSpace = std::static_pointer_cast<FunctionSpaceType>(this->data_.functionSpace());

  // merge local changes on the partitioned vector
  rightHandSide->setRepresentationGlobal();
  rightHandSide->startGhostManipulation();
  
  // get all entries
  std::vector<double> rhsValues;
  rightHandSide->getValuesWithGhosts(rhsValues);
  VLOG(1) << "extracted rhsValues (with ghosts): " << rhsValues;

  // initialize values to zero
  rightHandSide->zeroEntries();

  // also zero out the ghost buffer
  rightHandSide->zeroGhostBuffer();

  // setup arrays used for integration
  std::array<std::array<double,D>, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();
  EvaluationsArrayType evaluationsArray{};

  LOG(DEBUG) << "1D integration with " << QuadratureType::numberEvaluations() << " evaluations";
  LOG(DEBUG) << D << "D integration with " << QuadratureDD::numberEvaluations() << " evaluations";

  // set entries in rhs vector
  // loop over local elements
  for (element_no_t elementNo = 0; elementNo < functionSpace->nElementsLocal(); elementNo++)
  {
    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNosLocal = functionSpace->getElementDofNosLocal(elementNo);

    VLOG(2) << "element " << elementNo;

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
        // integrate value and set entry in stiffness matrix
        double integratedValue = integratedValues(i,j);

        double value = integratedValue * rhsValues[dofNosLocal[j]];
        VLOG(2) << "  dof pair (" << i<< "," <<j<< "), integrated value: " <<integratedValue << ", rhsValue[" << dofNosLocal[j]<< "]: " << rhsValues[dofNosLocal[j]] << " = " << value;

        rightHandSide->setValue(dofNosLocal[i], value, ADD_VALUES);
      }  // j
    }  // i
  }  // elementNo

  // merge local changes on the vector, parallel assembly
  rightHandSide->finishGhostManipulation();
 
}

};
