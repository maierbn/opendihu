#include "spatial_discretization/finite_element_method/04_rhs.h"

#include <Python.h>  // has to be the first included header
#include <array>
#include "quadrature/triangular_prism.h"

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
  typedef Quadrature::TriangularPrism<D,QuadratureType> QuadraturePrism;   // corner elements with triangles

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

  // if there is the extra vector that will hold rhs_active
  if (rightHandSideActive)
  {
    rightHandSideActive->zeroEntries();
    rightHandSideActive->setRepresentationGlobal();
    rightHandSideActive->startGhostManipulation();
  }

  // setup arrays used for integration
  std::array<std::array<double,D>, QuadratureDD::numberEvaluations()> samplingPointsHex = QuadratureDD::samplingPoints();
  std::array<std::array<double,D>, QuadraturePrism::numberEvaluations()> samplingPointsPrism = QuadraturePrism::samplingPoints();

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

    // determine if the element is a triangular prism at the corners or if it is a normal hex element
    ::Mesh::face_or_edge_t edge;
    bool isElementPrism = functionSpace->hasTriangleCorners()
      && functionSpace->meshPartition()->elementIsAtCorner(elementNoLocal, edge);

    // get number of sampling points
    int nSamplingPoints = samplingPointsHex.size();

    if (isElementPrism)
    {
      nSamplingPoints = samplingPointsPrism.size();
    }

    // compute integral
    for (unsigned int samplingPointIndex = 0; samplingPointIndex < nSamplingPoints; samplingPointIndex++)
    {
      // depending on the element type (hex or prism), get the sampling points
      std::array<double,D> xi;
      if (isElementPrism)
      {
        xi = samplingPointsPrism[samplingPointIndex];
      }
      else
      {
        xi = samplingPointsHex[samplingPointIndex];
      }
      // evaluate function to integrate at samplingPoints[i], write value to evaluations[i]

      // compute integration factor
      const std::array<Vec3,D> jacobian = functionSpace->computeJacobian(geometryValues, xi, elementNoLocal);
      double integrationFactor = MathUtility::computeIntegrationFactor(jacobian);

      Tensor2<D> inverseJacobian = functionSpace->getInverseJacobian(geometryValues, elementNoLocal, xi);

      std::array<VecD<D>,nDofsPerElement> gradPhiParameterSpace = functionSpace->getGradPhi(xi, elementNoLocal);

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
            double gradPhi = gradPhiWorldSpace[b];

            // active stress only contracts
            // If grad phi is positive, this means extension in the muscle tissue.
            // If ε>0, there is no value of the active stress. The virtual work -δW_ext = δW_int = σ_active : ε is therefore zero.
            // If ε<0, we have contraction. There exists a virtual work -δW_ext = δW_int = σ_active : ε.
            if (gradPhi < 0)
            {
              entryLa += gradPhi * sigma_ab;
            }
            VLOG(2) << "     " << samplingPointIndex << "," << dofIndexL << ", sigma_" << a << b << ": " << sigma_ab << ", grad phi: " << gradPhiWorldSpace[b];
          }
          evaluationsArray[samplingPointIndex][dofIndexL*nComponents + a] = entryLa * integrationFactor;
        }
      } // index L
    }  // function evaluations

    // integrate all values at once

    // depending on element type (triangular prism at the corners or normal hexahedral), use different quadrature
    EvaluationsType integratedValues;
    if (isElementPrism)
    {
      integratedValues = QuadraturePrism::computeIntegral(evaluationsArray);

      // set entries that are no real dofs in triangular prisms to zero
      const bool usesDofPairs = true;
      const bool isQuadraticElement0 = FunctionSpaceType::BasisFunction::getBasisOrder() == 2;
      const bool isQuadraticElement1 = FunctionSpaceType::BasisFunction::getBasisOrder() == 2;
      QuadraturePrism::adjustEntriesforPrism(integratedValues, edge, usesDofPairs,
                                             isQuadraticElement0, nComponents, isQuadraticElement1, nComponents);
    }
    else
    {
      integratedValues = QuadratureDD::computeIntegral(evaluationsArray);
    }

    VLOG(1) << "         evaluationsArray: " << evaluationsArray;
    VLOG(1) << "         integratedValues: " << integratedValues;

    // perform integration and add to entry in rhs vector
    for (int dofIndexL = 0; dofIndexL < nDofsPerElement; dofIndexL++)
    {
      int k = int(dofNosLocal[dofIndexL] / (this->data_.functionSpace()->meshPartition()->nNodesGlobal(0) * this->data_.functionSpace()->meshPartition()->nNodesGlobal(1)));
      VLOG(1) << "node k: " << k << ", f: " << integratedValues[dofIndexL*nComponents + 0] << " "
        << integratedValues[dofIndexL*nComponents + 1] << " " << integratedValues[dofIndexL*nComponents + 2]
        << ", sigma: " << activeStressValues[dofIndexL][MathUtility::sqr(nComponents)-1];

      // loop over components
      for (int componentNo = 0; componentNo < nComponents; componentNo++)
      {
        // there is no minus sign here
        double value = integratedValues[dofIndexL*nComponents + componentNo];

        VLOG(1) << "active rhs, dof " << dofNosLocal[dofIndexL] << ", component " << componentNo
          << ", add value " << value << " active stress: " << activeStressValues;
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
}

}  // namespace
