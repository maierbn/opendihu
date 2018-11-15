#include "postprocessing/streamline_tracer_base.h"

#include <algorithm>

namespace Postprocessing
{

template<typename FunctionSpace>
void StreamlineTracerBase<FunctionSpace>::
traceStreamline(Vec3 startingPoint, double direction, std::vector<Vec3> &points)
{
  const int D = FunctionSpace::dim();
  const int nDofsPerElement = FunctionSpace::nDofsPerElement();
  
  LOG(DEBUG) << "traceStreamline(startingPoint " << startingPoint << ", direction " << direction << ")";

  Vec3 currentPoint = startingPoint;
  element_no_t elementNo;
  int ghostMeshNo = -1;
  std::array<double,(unsigned long int)3> xi;
  
  // There are 2 implementations of streamline tracing.
  // The first one (useGradientField_) uses a precomputed gradient field that is interpolated linearly and the second uses the gradient directly from the Laplace solution_ field.
  // The first one seems more stable, because the gradient is zero and the position of the boundary conditions and should be used with a linear discretization of the potential field.
  // The second one is more accurate.

  std::array<Vec3,nDofsPerElement> elementalGradientValues;
  std::array<double,nDofsPerElement> elementalSolutionValues;
  std::array<Vec3,nDofsPerElement> geometryValues;

  bool startSearchInCurrentElement = false;

  // loop over length of streamline, avoid loops by limiting the number of iterations
  for (int iterationNo = 0; iterationNo <= maxNIterations_; iterationNo++)
  {
    if (iterationNo == maxNIterations_)
    {
      LOG(WARNING) << "streamline reached maximum number of iterations (" << maxNIterations_ << ")";
      points.clear();
      break;
    }

    if (iterationNo == 0)
    {
      startSearchInCurrentElement = false;
    }
    else
      startSearchInCurrentElement = true;

    bool positionFound = functionSpace_->findPosition(currentPoint, elementNo, ghostMeshNo, xi, startSearchInCurrentElement);

    // if no position was found, the streamline exits the domain
    if (!positionFound)
    {
      VLOG(2) << "streamline ends at iteration " << iterationNo << " because " << currentPoint << " is outside of domain";
      break;
    }

    // get values for element that are later needed to compute the gradient
    if (useGradientField_)
    {
      gradient_->getElementValues(elementNo, elementalGradientValues);
    }
    else
    {
      solution_->getElementValues(elementNo, elementalSolutionValues);

      // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
      functionSpace_->getElementGeometry(elementNo, geometryValues);
    }

    // get value of gradient
    Vec3 gradient;
    if (useGradientField_)
    {
      gradient = functionSpace_->template interpolateValueInElement<3>(elementalGradientValues, xi);
      VLOG(2) << "use gradient field";
    }
    else
    {
      // compute the gradient value in the current value
      Tensor2<D> inverseJacobian = functionSpace_->getInverseJacobian(geometryValues, elementNo, xi);
      gradient = functionSpace_->interpolateGradientInElement(elementalSolutionValues, inverseJacobian, xi);

      VLOG(2) << "use direct gradient";
    }

    // integrate streamline
    VLOG(2) << "  integrate from " << currentPoint << ", gradient: " << gradient << ", gradient normalized: " << MathUtility::normalized<3>(gradient)
      << ", lineStepWidth: " << lineStepWidth_;
    currentPoint = currentPoint + MathUtility::normalized<3>(gradient)*lineStepWidth_*direction;

    VLOG(2) << "              to " << currentPoint;

    points.push_back(currentPoint);
  }
}

};
