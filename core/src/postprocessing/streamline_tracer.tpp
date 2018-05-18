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
traceStreamlines()
{
  LOG(TRACE) << "traceStreamlines";
  
  // compute a gradient field from the solution
  problem_.data().solution().computeGradientField(data_.gradient());
  
  const int nDofsPerElement = DiscretizableInTimeType::BasisOnMesh::nDofsPerElement();
  std::vector<Vec3> nodePositions;
  
  // loop over seed points
  #pragma omp parallel for
  for (int seedPointNo = 0; seedPointNo != seedPositions_.size(); seedPointNo++)
  {
    // get starting point
    Vec3 currentPoint = seedPositions_[seedPointNo];
    nodePositions.push_back(currentPoint);
   
    // find out element 
    std::array<double,(unsigned long int)3> xi;
    element_no_t elementNo;
    std::array<Vec3,nDofsPerElement> elementalGradientValues;
    
    // find out initial element no and xi value
    bool positionFound = problem_.data().mesh()->findPosition(currentPoint, elementNo, xi);
    if (!positionFound)
    {
      LOG(ERROR) << "Seed point " << currentPoint << " is outside of domain.";
      continue;
    }
    
    // get gradient values for element 
    data_.gradient().getElementValues(elementNo, elementalGradientValues);
    
    for(;;)
    {
      
      // check if element_no is still valid 
      if (!problem_.data().mesh()->pointIsInElement(currentPoint, elementNo, xi))
      {
        bool positionFound = problem_.data().mesh()->findPosition(currentPoint, elementNo, xi);
        
        // if no position was found the streamline exists the domain
        if (!positionFound)
        {
          break;
        }
      }
      
      
      // get value of gradient 
      Vec3 gradient = problem_.data().mesh()->template interpolateValueInElement<3>(elementalGradientValues, xi);
      
      // integrate streamline
      currentPoint = currentPoint + gradient*lineStepWidth_;
      nodePositions.push_back(currentPoint);
    }
  }
}
 
};
