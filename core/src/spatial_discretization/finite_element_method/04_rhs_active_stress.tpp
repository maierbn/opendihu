#include "spatial_discretization/finite_element_method/04_rhs.h"

#include <Python.h>  // has to be the first included header
#include <array>

namespace SpatialDiscretization
{

template<typename FunctionSpaceType, typename QuadratureType, int nComponents>
void FiniteElementMethodRhs<FunctionSpaceType, QuadratureType, nComponents, Equation::Static::LinearElasticityActiveStress>::
setRightHandSide()
{
  LOG(TRACE) << "setRightHandSide LinearElasticityActiveStress";

  // compute rhs as f_La = -int sigma_ab * dphi_L(x)/dx_b dx

  const int D = FunctionSpaceType::dim();
  assert(D == nComponents);

  // define shortcuts for quadrature and basis
  typedef Quadrature::TensorProduct<D,QuadratureType> QuadratureDD;
  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();
  const int nUnknownsPerElement = nDofsPerElement*nComponents;
  typedef std::array<double,nUnknownsPerElement> EvaluationsType;
  typedef std::array<
            EvaluationsType,
            QuadratureDD::numberEvaluations()
          > EvaluationsArrayType;    // evaluations[nGP^D][nDofs][nDofs]

  // get shortcuts to variables
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> rightHandSide = this->data_.rightHandSide();
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents*nComponents>> activeStress = this->data_.activeStress();
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> rightHandSideActive = this->data_.rightHandSideActive();

  assert(activeStress != nullptr);

  std::shared_ptr<FunctionSpaceType> functionSpace = this->data_.functionSpace();

  // merge local changes on the partitioned vector
  rightHandSide->zeroEntries();
  rightHandSide->setRepresentationGlobal();
  rightHandSide->startGhostManipulation();

  if (rightHandSideActive)
  {
    rightHandSideActive->zeroEntries();
    rightHandSideActive->setRepresentationGlobal();
    rightHandSideActive->startGhostManipulation();
  }

  // setup arrays used for integration
  std::array<std::array<double,D>, QuadratureDD::numberEvaluations()> samplingPoints = QuadratureDD::samplingPoints();
  EvaluationsArrayType evaluationsArray{};

  // set entries in rhs vector
  // loop over local elements
  for (element_no_t elementNoLocal = 0; elementNoLocal < functionSpace->nElementsLocal(); elementNoLocal++)
  {
    // get indices of element-local dofs
    std::array<dof_no_t,nDofsPerElement> dofNosLocal = functionSpace->getElementDofNosLocal(elementNoLocal);

    // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
    std::array<Vec3,nDofsPerElement> geometryValues;
    functionSpace->getElementGeometry(elementNoLocal, geometryValues);

    std::array<VecD<nComponents*nComponents>,nDofsPerElement> activeStressValues;
    activeStress->getElementValues(elementNoLocal, activeStressValues);

    // compute integral
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < samplingPoints.size(); samplingPointIndex++)
    {
      // evaluate function to integrate at samplingPoints[i], write value to evaluations[i]
      std::array<double,D> xi = samplingPoints[samplingPointIndex];

      // compute integration factor
      const std::array<Vec3,D> jacobian = FunctionSpaceType::computeJacobian(geometryValues, xi);
      double integrationFactor = MathUtility::computeIntegrationFactor(jacobian);

      Tensor2<D> inverseJacobian = functionSpace->getInverseJacobian(geometryValues, elementNoLocal, xi);

      std::array<VecD<D>,nDofsPerElement> gradPhiParameterSpace = functionSpace->getGradPhi(xi);

      // loop over dofs of element (index L in formula)
      for (int dofIndexL = 0; dofIndexL < nDofsPerElement; dofIndexL++)
      {
        // compute grad phi with respect to world space
        VecD<D> gradPhiWorldSpace = inverseJacobian * gradPhiParameterSpace[dofIndexL];

        // compute entry f_active_La = int_Omega σ_ab ∂phi_L/∂x_b
        for (int a = 0; a < D; a++)
        {
          double entryLa = 0;
          for (int b = 0; b < D; b++)
          {
            double sigma_ab = activeStressValues[dofIndexL][a*D + b];
            entryLa += gradPhiWorldSpace[b] * sigma_ab;
          }
          evaluationsArray[samplingPointIndex][dofIndexL*nComponents + a] = entryLa * integrationFactor;
        }
      } // index L
    }  // function evaluations

    // integrate all values at once
    EvaluationsType integratedValues = QuadratureDD::computeIntegral(evaluationsArray);

    // perform integration and add to entry in rhs vector
    for (int dofIndexL = 0; dofIndexL < nDofsPerElement; dofIndexL++)
    {
      // loop over components
      for (int componentNo = 0; componentNo < nComponents; componentNo++)
      {
        double value = -integratedValues[dofIndexL*nComponents + componentNo];

        LOG(DEBUG) << "rhs, dof " << dofNosLocal[dofIndexL] << ", component " << componentNo << ", add value " << value;
        rightHandSide->setValue(componentNo, dofNosLocal[dofIndexL], value, ADD_VALUES);
        if(rightHandSideActive)
          rightHandSideActive->setValue(componentNo, dofNosLocal[dofIndexL], value, ADD_VALUES);
      }
    } // index L
  }  // elementNoLocal

  // merge local changes on the vector, parallel assembly
  rightHandSide->finishGhostManipulation();
  if (rightHandSideActive)
    rightHandSideActive->finishGhostManipulation();
  LOG(DEBUG) << "set rhs: " << *rightHandSide;

}

}  // namespace
