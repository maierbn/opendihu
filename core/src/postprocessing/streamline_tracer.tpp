#include "postprocessing/streamline_tracer.h"

#include <petscvec.h>
#include "utility/python_utility.h"

namespace Postprocessing
{

template<typename DiscretizableInTimeType>
StreamlineTracer<DiscretizableInTimeType>::
StreamlineTracer(DihuContext context) :
  context_(context["StreamlineTracer"]), problem_(context_), data_(context_)
{
  LOG(TRACE) << "StreamlineTracer::StreamlineTracer()";

  specificSettings_ = context_.getPythonConfig();
  VLOG(2) << "in StreamlineTracer(), specificSettings_: " << specificSettings_;
  outputWriterManager_.initialize(specificSettings_);

  lineStepWidth_ = PythonUtility::getOptionDouble(specificSettings_, "lineStepWidth", 1e-2);
  maxNIterations_ = PythonUtility::getOptionInt(specificSettings_, "maxIterations", 100000, PythonUtility::Positive);
  useGradientField_ = PythonUtility::getOptionBool(specificSettings_, "useGradientField_", false);
   
  // get the first seed position from the list
  PyObject *pySeedPositions = PythonUtility::getOptionListBegin<PyObject *>(specificSettings_, "seedPoints");

  // loop over other entries of list
  for (;
      !PythonUtility::getOptionListEnd(specificSettings_, "seedPoints");
      PythonUtility::getOptionListNext<PyObject *>(specificSettings_, "seedPoints", pySeedPositions))
  {
    Vec3 seedPosition = PythonUtility::convertFromPython<Vec3>(pySeedPositions);
    seedPositions_.push_back(seedPosition);
  }
}

template<typename DiscretizableInTimeType>
void StreamlineTracer<DiscretizableInTimeType>::
initialize()
{
  LOG(TRACE) << "StreamlineTracer::initialize";

  // initialize the problem
  problem_.initialize();

  // initialize streamline tracer data object
  data_.setBaseData(std::make_shared<typename DiscretizableInTimeType::Data>(problem_.data()));
  data_.initialize();
}

template<typename DiscretizableInTimeType>
void StreamlineTracer<DiscretizableInTimeType>::
run()
{
  initialize();

  // call the method of the underlying problem
  problem_.run();

  // do the tracing
  traceStreamlines();

  // output
  outputWriterManager_.writeOutput(data_);
}

template<typename DiscretizableInTimeType>
void StreamlineTracer<DiscretizableInTimeType>::
traceStreamline(element_no_t initialElementNo, std::array<double,(unsigned long int)3> xi, Vec3 startingPoint, double direction, std::vector<Vec3> &points)
{
  const int D = DiscretizableInTimeType::BasisOnMesh::dim();
  
  const int nDofsPerElement = DiscretizableInTimeType::BasisOnMesh::nDofsPerElement();
  
  Vec3 currentPoint = startingPoint;
  element_no_t elementNo = initialElementNo;
  
    std::array<Vec3,nDofsPerElement> elementalGradientValues;
    std::array<double,nDofsPerElement> elementalSolutionValues;
    std::array<Vec3,nDofsPerElement> geometryValues;

   // get gradient values for element
   // There are 2 implementations of streamline tracing. 
   // The first one (useGradientField_) uses a precomputed gradient field that is interpolated linearly and the second uses the gradient directly from the Laplace solution field. 
   // The first one seems more stable, because the gradient is zero and the position of the boundary conditions and should be used with a linear discretization of the potential field. 
   // The second one is more accurate.
   if (useGradientField_)
   {
     // use the precomputed gradient field
     data_.gradient().getElementValues(elementNo, elementalGradientValues);
   }
   else 
   {
     // get the local gradient value at the current position
     problem_.data().solution().getElementValues(elementNo, elementalSolutionValues);

     // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
     problem_.data().mesh()->getElementGeometry(elementNo, geometryValues);
   }
   
   VLOG(2) << "streamline starts in element " << elementNo;
   
   // loop over length of streamline, avoid loops by limiting the number of iterations
   for(int iterationNo = 0; iterationNo <= maxNIterations_; iterationNo++)
   {
     if (iterationNo == maxNIterations_)
     {
       LOG(WARNING) << "streamline reached maximum number of iterations (" << maxNIterations_ << ")";
       points.clear();
       break;
     }
    
     // check if element_no is still valid
     if (!problem_.data().mesh()->pointIsInElement(currentPoint, elementNo, xi))
     {
       bool positionFound = problem_.data().mesh()->findPosition(currentPoint, elementNo, xi);

       // if no position was found, the streamline exits the domain
       if (!positionFound)
       {
         VLOG(2) << "streamline ends at iteration " << iterationNo << " because " << currentPoint << " is outside of domain";
         break;
       }
           
       // get values for element that are later needed to compute the gradient
       if (useGradientField_)      
       {
         data_.gradient().getElementValues(elementNo, elementalGradientValues);
       }
       else 
       {
         problem_.data().solution().getElementValues(elementNo, elementalSolutionValues);
           
         // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
         problem_.data().mesh()->getElementGeometry(elementNo, geometryValues);
       }
       
       VLOG(2) << "streamline enters element " << elementNo;
     }

     // get value of gradient
     Vec3 gradient;
     if (useGradientField_)
     {       
       gradient = problem_.data().mesh()->template interpolateValueInElement<3>(elementalGradientValues, xi);
       VLOG(2) << "use gradient field";
     }
     else 
     {
       // compute the gradient value in the current value
       Tensor2<D> inverseJacobian = problem_.data().mesh()->getInverseJacobian(geometryValues, elementNo, xi);
       gradient = problem_.data().mesh()->interpolateGradientInElement(elementalSolutionValues, inverseJacobian, xi);
       
       VLOG(2) << "use direct gradient";
     }
     
     // integrate streamline
     VLOG(2) << "  integrate from " << currentPoint << ", gradient: " << gradient << ", gradient normalized: " << MathUtility::normalized<3>(gradient) << ", lineStepWidth: " << lineStepWidth_;
     currentPoint = currentPoint + MathUtility::normalized<3>(gradient)*lineStepWidth_*direction;
     
     VLOG(2) << "              to " << currentPoint;
     
     points.push_back(currentPoint);
   }
}


template<typename DiscretizableInTimeType>
void StreamlineTracer<DiscretizableInTimeType>::
traceStreamlines()
{
  LOG(TRACE) << "traceStreamlines";

  // compute a gradient field from the solution
  problem_.data().solution().computeGradientField(data_.gradient());
 
  std::array<double,(unsigned long int)3> xi;
  std::vector<std::vector<Vec3>> nodePositions(seedPositions_.size());
  const int nSeedPoints = seedPositions_.size();

  LOG(DEBUG) << "trace streamline, seedPositions: " << seedPositions_;
  
  // loop over seed points
  //#pragma omp parallel for shared(nodePositions)
  for (int seedPointNo = 0; seedPointNo < nSeedPoints; seedPointNo++)
  {
    LOG(DEBUG) << " seed point " << seedPointNo;
   
    // get starting point
    Vec3 startingPoint = seedPositions_[seedPointNo];
    nodePositions[seedPointNo].push_back(startingPoint);

    element_no_t initialElementNo = 0;
    
    // find out initial element no and xi value where the current Point lies
    bool positionFound = problem_.data().mesh()->findPosition(startingPoint, initialElementNo, xi);
    if (!positionFound)
    {
      LOG(ERROR) << "Seed point " << startingPoint << " is outside of domain.";
      continue;
    }
    
    // trace streamline forwards
    std::vector<Vec3> forwardPoints;
    traceStreamline(initialElementNo, xi, startingPoint, 1.0, forwardPoints);
    
    // trace streamline backwards
    std::vector<Vec3> backwardPoints;
    traceStreamline(initialElementNo, xi, startingPoint, -1.0, backwardPoints);
  
    // copy collected points to result vector 
    nodePositions[seedPointNo].insert(nodePositions[seedPointNo].begin(), backwardPoints.rbegin(), backwardPoints.rend());
    nodePositions[seedPointNo].insert(nodePositions[seedPointNo].end(), startingPoint);
    nodePositions[seedPointNo].insert(nodePositions[seedPointNo].end(), forwardPoints.begin(), forwardPoints.end());
  }
  
  // create 1D meshes of streamline from collected node positions
  std::ofstream file("streamlines.csv");
  if (!file.is_open())
    LOG(WARNING) << "Could not open streamlines.csv for writing";
  
  for (int seedPointNo = 0; seedPointNo != seedPositions_.size(); seedPointNo++)
  {
    for (std::vector<Vec3>::const_iterator iter = nodePositions[seedPointNo].begin(); iter != nodePositions[seedPointNo].end(); iter++)
    {
      Vec3 point = *iter;
      file << point[0] << ";" << point[1] << ";" << point[2] << ";";
    }
    file << "\n";
  }
  file.close();
  
  // create new meshes, one for each streamline 
  for (int seedPointNo = 0; seedPointNo != seedPositions_.size(); seedPointNo++)
  {
    LOG(DEBUG) << "seed point " << seedPointNo << ", number node positions: " << nodePositions[seedPointNo].size();
    this->data_.createFibreMesh(nodePositions[seedPointNo]);
  }
}

};
