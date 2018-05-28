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

#define USE_GRADIENT_FIELD
  
#ifndef USE_GRADIENT_FIELD
  const int D = DiscretizableInTimeType::BasisOnMesh::dim();
#endif
  
  const int nDofsPerElement = DiscretizableInTimeType::BasisOnMesh::nDofsPerElement();
  std::vector<std::vector<Vec3>> nodePositions(seedPositions_.size());
  const int nSeedPoints = seedPositions_.size();

  LOG(DEBUG) << "trace streamline, seedPositions: " << seedPositions_;
  
  // loop over seed points
  //#pragma omp parallel for shared(nodePositions)
  for (int seedPointNo = 0; seedPointNo < nSeedPoints; seedPointNo++)
  {
    LOG(DEBUG) << " seed point " << seedPointNo;
   
    // get starting point
    Vec3 currentPoint = seedPositions_[seedPointNo];
    nodePositions[seedPointNo].push_back(currentPoint);

    std::array<double,(unsigned long int)3> xi;
    element_no_t elementNo = 0;
    
#ifdef USE_GRADIENT_FIELD    
    std::array<Vec3,nDofsPerElement> elementalGradientValues;
#else    
    std::array<double,nDofsPerElement> elementalSolutionValues;
    std::array<Vec3,nDofsPerElement> geometryValues;
#endif

    // find out initial element no and xi value where the current Point lies
    bool positionFound = problem_.data().mesh()->findPosition(currentPoint, elementNo, xi);
    if (!positionFound)
    {
      LOG(ERROR) << "Seed point " << currentPoint << " is outside of domain.";
      continue;
    }

    // get gradient values for element
#ifdef USE_GRADIENT_FIELD    
    data_.gradient().getElementValues(elementNo, elementalGradientValues);
#else    
    problem_.data().solution().getElementValues(elementNo, elementalSolutionValues);

    // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
    problem_.data().mesh()->getElementGeometry(elementNo, geometryValues);
#endif
    
    VLOG(2) << "streamline starts in element " << elementNo;
        
    // loop over length of streamline
    for(int iterationNo = 0; iterationNo <= 1000000; iterationNo++)
    {
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
            
        // get gradient values for element
#ifdef USE_GRADIENT_FIELD          
        data_.gradient().getElementValues(elementNo, elementalGradientValues);
#else
        problem_.data().solution().getElementValues(elementNo, elementalSolutionValues);
            
        // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
        problem_.data().mesh()->getElementGeometry(elementNo, geometryValues);
#endif
        
        VLOG(2) << "streamline enters element " << elementNo;
      }

      // get value of gradient

#ifdef USE_GRADIENT_FIELD          
      Vec3 gradient = problem_.data().mesh()->template interpolateValueInElement<3>(elementalGradientValues, xi);
#else
      Tensor2<D> inverseJacobian = problem_.data().mesh()->getInverseJacobian(geometryValues, elementNo, xi);
      Vec3 gradient = problem_.data().mesh()->interpolateGradientInElement(elementalSolutionValues, inverseJacobian, xi);
#endif
      // integrate streamline
      
      VLOG(2) << "  integrate from " << currentPoint << ", gradient: " << gradient << ", gradient normalized: " << MathUtility::normalized<3>(gradient) << ", lineStepWidth: " << lineStepWidth_;
      currentPoint = currentPoint + MathUtility::normalized<3>(gradient)*lineStepWidth_;
      
      VLOG(2) << "              to " << currentPoint;
      
      nodePositions[seedPointNo].push_back(currentPoint);
      
    }
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
    this->data_.createFibreMesh(nodePositions[seedPointNo]);
  }
}

};
