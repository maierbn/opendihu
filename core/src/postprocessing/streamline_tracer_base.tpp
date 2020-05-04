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
  
  LOG(DEBUG) << "traceStreamline(startingPoint " << startingPoint << ", direction " << direction << ", maxNIterations_: " << maxNIterations_ << "), useGradientField_: " << useGradientField_;

  //functionSpace_->debugOutputGhostMeshSet();

  Vec3 currentPoint = startingPoint;
  element_no_t elementNo = 0;
  int ghostMeshNo = -1;
  std::array<double,(unsigned long int)3> xi;
  std::shared_ptr<FunctionSpace> functionSpace = functionSpace_;  //< the function space to use, this can be set to one of the ghost meshes
  
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
    {
      startSearchInCurrentElement = true;
    }

    VLOG(1) << "startSearchInCurrentElement=" << startSearchInCurrentElement << ", elementNo: " << elementNo << " ghostMeshNo: " << ghostMeshNo;
    double residual;
    bool searchedAllElements;

    // look for the element and xi value of the currentPoint, also considers ghost meshes if they are set
    bool positionFound = functionSpace_->findPosition(currentPoint, elementNo, ghostMeshNo, xi, startSearchInCurrentElement, residual, searchedAllElements);

    // if no position was found, the streamline exits the domain
    if (!positionFound)
    {
      VLOG(2) << "streamline ends at iteration " << iterationNo << " because " << currentPoint << " is outside of domain";
      break;
    }

    VLOG(1) << " findPosition returned ghostMeshNo " << ghostMeshNo << ", elementNo " << elementNo << ", xi " << xi;

    // get values for element that are later needed to compute the gradient

    // if the streamline passes a normal element
    if (ghostMeshNo == -1)
    {
      VLOG(1) << "use normal mesh";

      functionSpace = functionSpace_;

      if (useGradientField_)
      {
        gradient_->getElementValues(elementNo, elementalGradientValues);
      }
      else
      {
        solution_->getElementValues(elementNo, elementalSolutionValues);

        // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
        functionSpace->getElementGeometry(elementNo, geometryValues);
      }
    }
    else    // if the streamline is in an element of a ghost mesh
    {
      VLOG(1) << "use ghost mesh";

      // use ghost mesh as current function space
      functionSpace = ghostMesh_[ghostMeshNo];

      if (useGradientField_)
      {
        ghostMeshGradient_[ghostMeshNo]->getElementValues(elementNo, elementalGradientValues);
      }
      else
      {
        ghostMeshSolution_[ghostMeshNo]->getElementValues(elementNo, elementalSolutionValues);

        // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
        functionSpace->getElementGeometry(elementNo, geometryValues);
      }
    }

    // get value of gradient
    Vec3 gradient;
    if (useGradientField_)
    {
      gradient = functionSpace->template interpolateValueInElement<3>(elementalGradientValues, xi);
      VLOG(2) << "use gradient field";
    }
    else
    {
      // compute the gradient value in the current value
      Tensor2<D> inverseJacobian = functionSpace->getInverseJacobian(geometryValues, elementNo, xi);
      gradient = functionSpace->interpolateGradientInElement(elementalSolutionValues, inverseJacobian, xi);

      VLOG(2) << "use direct gradient";
    }

    if (fabs(gradient[0] + gradient[1] + gradient[2]) < 1e-15)
    {
      LOG(ERROR) << "Gradient at element " << elementNo << ", xi " << xi << " is zero!";
      if (!useGradientField_)
      {
        Tensor2<D> inverseJacobian = functionSpace->getInverseJacobian(geometryValues, elementNo, xi);
        LOG(DEBUG) << "geometryValues: " << geometryValues << ", inverseJacobian: " << inverseJacobian << ", elementalSolutionValues: " << elementalSolutionValues;
      }
      else
      {
        LOG(DEBUG) << "elementalGradientValues: " << elementalGradientValues;
      }
      MPI_Abort(functionSpace_->meshPartition()->mpiCommunicator(), -1);
    }

    // integrate streamline
    VLOG(2) << "  integrate from " << currentPoint << ", gradient: " << gradient << ", gradient normalized: " << MathUtility::normalized<3>(gradient)
      << ", lineStepWidth: " << lineStepWidth_;
    currentPoint = currentPoint + MathUtility::normalized<3>(gradient)*lineStepWidth_*direction;

    VLOG(2) << "              to " << currentPoint;

    points.push_back(currentPoint);
  }

  if (points.empty())
  {
    LOG(DEBUG) << "traced streamline is completely empty";
  }
  else if (points.size() == 1)
  {
    LOG(DEBUG) << "traced streamline has 1 point: " << points[0];
  }
  else
  {
    LOG(DEBUG) << "traced streamline has " << points.size() << " points, start: " << points[0] << ", end: " << points[points.size()-1];
  }
}

}  // namespace
