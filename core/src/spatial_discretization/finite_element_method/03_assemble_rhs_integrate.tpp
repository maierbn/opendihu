#include "spatial_discretization/finite_element_method/03_assemble_rhs.h"

#include <Python.h>
#include <memory>
#include <vector>
#include <petscsys.h>
#include <array>
#include <functional>

#include "quadrature/tensor_product.h"
#include "function_space/function_space.h"
#include "spatial_discretization/finite_element_method/integrand/integrand_mass_matrix.h"
#include "field_variable/field_variable.h"
#include "utility/vector_operators.h"


namespace SpatialDiscretization
{

// 1D,2D,3D rhs vector of Deformable mesh
template<typename FunctionSpaceType, typename QuadratureType, int nComponents, typename Term, typename Dummy>
void AssembleRightHandSide<FunctionSpaceType, QuadratureType, nComponents, Term, Dummy>::
multiplyRightHandSideWithMassMatrix()
{
  const int D = FunctionSpaceType::dim();
  LOG(TRACE) << "multiplyRightHandSideWithMassMatrix " << D << "D";

  // define shortcuts for integrator and basis
  typedef Quadrature::TensorProduct<D,QuadratureType> QuadratureDD;
  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();
  const int nUnknownsPerElement = nDofsPerElement*nComponents;
  typedef MathUtility::Matrix<nUnknownsPerElement,nUnknownsPerElement,double_v_t> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureDD::numberEvaluations()
          > EvaluationsArrayType;    // evaluations[nGP^D][nDofs][nDofs]

  // setup arrays used for integration
  std::array<std::array<double,D>, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();
  EvaluationsArrayType evaluationsArray{};

  LOG(DEBUG) << "1D integration with " << QuadratureType::numberEvaluations() << " evaluations";
  LOG(DEBUG) << D << "D integration with " << QuadratureDD::numberEvaluations() << " evaluations";

  // initialize variables
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> rightHandSide = this->data_.rightHandSide();

  std::shared_ptr<FunctionSpaceType> functionSpace = std::static_pointer_cast<FunctionSpaceType>(this->data_.functionSpace());

  // merge local changes on the partitioned vector
  rightHandSide->setRepresentationGlobal();
  rightHandSide->startGhostManipulation();
  
  // get all entries
  std::vector<VecD<nComponents>> rhsValues;
  rightHandSide->getValuesWithGhosts(rhsValues);
  VLOG(1) << "extracted rhsValues (with ghosts): " << rhsValues;

  // initialize values to zero
  rightHandSide->zeroEntries();

  // also zero out the ghost buffer
  rightHandSide->zeroGhostBuffer();

  // set entries in rhs vector

  element_no_t nElementsLocal = functionSpace->nElementsLocal();

  // initialize values to zero
  // loop over elements, always 4 elements at once using the vectorized functions
  for (int elementNoLocal = 0; elementNoLocal < nElementsLocal; elementNoLocal += nVcComponents)
  {

#ifdef USE_VECTORIZED_FE_MATRIX_ASSEMBLY
    // get indices of elementNos that should be handled in the current iterations,
    // this is, e.g.
    //    [10,11,12,13,-1,-1,-1,-1] (if nVcComponents==4 and nElementsLocal > 13)
    // or [10,11,12,-1,-1,-1,-1,-1] (if nVcComponents==4 and nElementsLocal == 13)

    dof_no_v_t elementNoLocalv([elementNoLocal, nElementsLocal](int i)
    {
      return (i >= nVcComponents || elementNoLocal+i >= nElementsLocal? -1: elementNoLocal+i);
    });

    // here, elementNoLocalv is the list of indices of the current iteration, e.g. [10,11,12,13,-1,-1,-1,-1]
    // elementNoLocal is the first entry of elementNoLocalv
#else
    int elementNoLocalv = elementNoLocal;
#endif

    // get indices of element-local dofs
    std::array<dof_no_v_t,nDofsPerElement> dofNosLocal = functionSpace->getElementDofNosLocal(elementNoLocalv);

    VLOG(2) << "element " << elementNoLocalv;

    // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
    std::array<Vec3_v_t,FunctionSpaceType::nDofsPerElement()> geometry;
    functionSpace->getElementGeometry(elementNoLocalv, geometry);

    // compute integral
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // evaluate function to integrate at samplingPoints[i*2], write value to evaluations[i]
      std::array<double,D> xi = samplingPoints[samplingPointIndex];

      // compute the 3xD jacobian of the parameter space to world space mapping
      auto jacobian = FunctionSpaceType::computeJacobian(geometry, xi);

      // get evaluations of integrand which is defined in another class
      evaluationsArray[samplingPointIndex] = IntegrandMassMatrix<D,EvaluationsType,FunctionSpaceType,nComponents,double_v_t,dof_no_v_t,Term>::
        evaluateIntegrand(jacobian,xi);

    }  // function evaluations

    // integrate all values for the (i,j) dof pairs at once
    EvaluationsType integratedValues = QuadratureDD::computeIntegral(evaluationsArray);

    // perform integration and add to entry in rhs vector
    for (int i = 0; i < nDofsPerElement; i++)
    {
      for (int j = 0; j < nDofsPerElement; j++)
      {
        // loop over components (1,...,D for solid mechanics)
        for (int rowComponentNo = 0; rowComponentNo < nComponents; rowComponentNo++)
        {
          for (int columnComponentNo = 0; columnComponentNo < nComponents; columnComponentNo++)
          {
            // integrate value and set entry in stiffness matrix
            double_v_t integratedValue = integratedValues(i*nComponents + rowComponentNo, j*nComponents + columnComponentNo);

            // getValuesAtIndices replaces operator[] (i.e. "rhsValues[dofNosLocal[j]]") and is necessary for it to work also with Vc::double_v
            double_v_t value = integratedValue * getValuesAtIndices(rhsValues,dofNosLocal[j])[columnComponentNo];
            VLOG(2) << "  dof pair (" << i<< "," <<j<< "), integrated value: " <<integratedValue << ", rhsValue[" << dofNosLocal[j]<< "]: " << getValuesAtIndices(rhsValues,dofNosLocal[j]) << " = " << value;

            rightHandSide->setValue(rowComponentNo, dofNosLocal[i], value, ADD_VALUES);
          }
        }
      }  // j
    }  // i
  }  // elementNoLocalv

  // merge local changes on the vector, parallel assembly
  rightHandSide->finishGhostManipulation();
 
}

}  // namespace
