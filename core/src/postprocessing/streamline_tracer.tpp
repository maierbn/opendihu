#include "postprocessing/streamline_tracer.h"

#include <algorithm>
#include <petscvec.h>

#include "utility/python_utility.h"

namespace Postprocessing
{

template<typename DiscretizableInTimeType>
StreamlineTracer<DiscretizableInTimeType>::
StreamlineTracer(DihuContext context) :
  context_(context["StreamlineTracer"]), problem_(context_), data_(context_), specificSettings_(context_.getPythonConfig())
{
  LOG(TRACE) << "StreamlineTracer::StreamlineTracer()";

  VLOG(2) << "in StreamlineTracer(), specificSettings_: " << specificSettings_;
  outputWriterManager_.initialize(context_, specificSettings_);

  this->lineStepWidth_ = specificSettings_.getOptionDouble("lineStepWidth", 1e-2, PythonUtility::Positive);
  this->maxNIterations_ = specificSettings_.getOptionInt("maxIterations", 100000, PythonUtility::Positive);
  this->useGradientField_ = specificSettings_.getOptionBool("useGradientField", false);

  targetElementLength_ = specificSettings_.getOptionDouble("targetElementLength", 0.0, PythonUtility::Positive);
  targetLength_ = specificSettings_.getOptionDouble("targetLength", 0.0, PythonUtility::Positive);
  discardRelativeLength_ = specificSettings_.getOptionDouble("discardRelativeLength", 0.0, PythonUtility::Positive);
  csvFilename_ = specificSettings_.getOptionString("csvFilename", "");
  csvFilenameBeforePostprocessing_ = specificSettings_.getOptionString("csvFilenameBeforePostprocessing", "");
  
  // get the first seed position from the list
  PyObject *pySeedPositions = specificSettings_.getOptionListBegin<PyObject *>("seedPoints");

  // loop over other entries of list
  for (;
      !specificSettings_.getOptionListEnd("seedPoints");
      specificSettings_.getOptionListNext<PyObject *>("seedPoints", pySeedPositions))
  {
    Vec3 seedPosition = PythonUtility::convertFromPython<Vec3>::get(pySeedPositions);
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

  // initialize values in base class
  this->functionSpace_ = problem_.data().functionSpace();
  this->solution_ = problem_.data().solution();
  this->gradient_ = data_.gradient();
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
  problem_.data().solution()->computeGradientField(data_.gradient());
 
  const int nSeedPoints = seedPositions_.size();
  //const int nSeedPoints = 1;
  std::vector<std::vector<Vec3>> streamlines(nSeedPoints);

  LOG(DEBUG) << "trace streamline, seedPositions: " << seedPositions_;
  
  // loop over seed points
  //#pragma omp parallel for shared(streamlines)
  for (int seedPointNo = 0; seedPointNo < nSeedPoints; seedPointNo++)
  {
    // get starting point
    Vec3 startingPoint = seedPositions_[seedPointNo];
    
    // trace streamline forwards
    std::vector<Vec3> forwardPoints;
    this->traceStreamline(startingPoint, 1.0, forwardPoints);
    
    if (forwardPoints.empty())  // if there was not even the first point found
    {
      LOG(ERROR) << "Seed point " << startingPoint << " is outside of domain.";
      continue;
    }

    // trace streamline backwards
    std::vector<Vec3> backwardPoints;
    this->traceStreamline(startingPoint, -1.0, backwardPoints);
  
    // copy collected points to result vector, note avoiding this additional copy-step is not really possible, since it would require a push_front which is only efficient with lists, but we need a vector here
    streamlines[seedPointNo].insert(streamlines[seedPointNo].begin(), backwardPoints.rbegin(), backwardPoints.rend());
    streamlines[seedPointNo].insert(streamlines[seedPointNo].end(), startingPoint);
    streamlines[seedPointNo].insert(streamlines[seedPointNo].end(), forwardPoints.begin(), forwardPoints.end());
    
    LOG(DEBUG) << " seed point " << seedPointNo << ", " << streamlines[seedPointNo].size() << " points";
  }
  
  // create 1D meshes of streamline from collected node positions
  if (!csvFilenameBeforePostprocessing_.empty())
  {
    std::ofstream file(csvFilenameBeforePostprocessing_, std::ios::out | std::ios::binary | std::ios::trunc);
    if (!file.is_open())
      LOG(WARNING) << "Could not open \"" << csvFilenameBeforePostprocessing_ << "\" for writing";
    
    for (int streamlineNo = 0; streamlineNo != streamlines.size(); streamlineNo++)
    {
      for (std::vector<Vec3>::const_iterator iter = streamlines[streamlineNo].begin(); iter != streamlines[streamlineNo].end(); iter++)
      {
        Vec3 point = *iter;
        file << point[0] << ";" << point[1] << ";" << point[2] << ";";
      }
      file << "\n";
    }
    file.close();
    LOG(INFO) << "File \"" << csvFilenameBeforePostprocessing_ << "\" written.";
  }
  
  // coarsen streamlines and drop too small streamlines
  postprocessStreamlines(streamlines);
 
  LOG(DEBUG) << "number streamlines after postprocessStreamlines: " << streamlines.size();
  
  // create 1D meshes of streamline from collected node positions
  if (!csvFilename_.empty())
  {
    std::ofstream file(csvFilename_, std::ios::out | std::ios::binary | std::ios::trunc);
    if (!file.is_open())
      LOG(WARNING) << "Could not open \"" << csvFilename_ << "\" for writing";
    
    for (int streamlineNo = 0; streamlineNo != streamlines.size(); streamlineNo++)
    {
      for (std::vector<Vec3>::const_iterator iter = streamlines[streamlineNo].begin(); iter != streamlines[streamlineNo].end(); iter++)
      {
        Vec3 point = *iter;
        file << point[0] << ";" << point[1] << ";" << point[2] << ";";
      }
      file << "\n";
    }
    file.close();
    LOG(INFO) << "File \"" << csvFilename_ << "\" written.";
  }
  
  // create new meshes, one for each streamline 
  for (int streamlineNo = 0; streamlineNo != streamlines.size(); streamlineNo++)
  {
    LOG(DEBUG) << "seed point " << streamlineNo << ", number node positions: " << streamlines[streamlineNo].size();
    this->data_.createfiberMesh(streamlines[streamlineNo]);
  }
}

template<typename DiscretizableInTimeType>
void StreamlineTracer<DiscretizableInTimeType>::
postprocessStreamlines(std::vector<std::vector<Vec3>> &streamlines)
{
  std::vector<double> lengths(streamlines.size());
 
  // compute length of each streamline 
  int i = 0;
  // loop over streamlines
  for (std::vector<std::vector<Vec3>>::iterator streamlinesIter = streamlines.begin(); streamlinesIter != streamlines.end(); streamlinesIter++, i++)
  {
    lengths[i] = 0.0;
    
    Vec3 lastPoint;
    bool firstPoint = true;
    int pointNo = 0;
    
    // loop over points of streamline
    for (std::vector<Vec3>::iterator pointsIter = streamlinesIter->begin(); pointsIter != streamlinesIter->end(); pointsIter++, pointNo++)
    {
      if (!firstPoint)
      {
        lengths[i] += MathUtility::distance<3>(*pointsIter, lastPoint);
      }
      firstPoint = false;
      lastPoint = *pointsIter;
    }
  }
  
  LOG(DEBUG) << " lengths of streamlines: " << lengths;
  
  // sort length
  std::vector<double> lengthsSorted(lengths);
  std::sort(lengthsSorted.begin(), lengthsSorted.end());
  
  // get median 
  double medianLength = lengthsSorted[lengthsSorted.size()/2];
  double maximumLength = lengthsSorted[lengthsSorted.size()-1];
  LOG(INFO) << "The median length of the streamlines is " << medianLength 
    << ", the maximum length of the " << lengthsSorted.size() << " streamlines is " << maximumLength << ".";
    
  if (discardRelativeLength_ != 0.0)
  {
   
    // clear streamlines that are shorter than discardRelativeLength_
    // loop over streamlines
    i = 0;
    for (std::vector<std::vector<Vec3>>::iterator streamlinesIter = streamlines.begin(); streamlinesIter != streamlines.end(); streamlinesIter++, i++)
    { 
      if (lengths[i] < discardRelativeLength_*medianLength) 
      {
        LOG(INFO) << "Discarding streamline no. " << i << " with length " << lengths[i] << " (Threshold " << discardRelativeLength_*medianLength << ").";
        streamlinesIter->clear();
      }
    }
    
    auto lastValidStreamline = std::remove_if (streamlines.begin(), streamlines.end(), 
                                              [](const std::vector<Vec3> &a)-> bool{return a.empty();});
    
    // remove previously cleared streamlines
    streamlines.erase(
       lastValidStreamline, 
       streamlines.end()
    );
  }
  
  // compute scale factor that scales streamlines to targetLength
  double scalingFactor = 1.0;
  if (targetLength_ != 0)
  {
    scalingFactor = targetLength_/maximumLength;
    LOG(INFO) << "Scaling factor for streamlines: " << scalingFactor;
  }
  
  // resample streamlines
  if (targetElementLength_ != 0.0 && targetElementLength_ != this->lineStepWidth_)
  {
    // loop over streamlines
    for (int i = 0; i < streamlines.size(); i++)
    {
      std::vector<Vec3> &currentStreamline = streamlines[i];
      
      if (currentStreamline.empty())
      {
        LOG(DEBUG) << "Streamline is empty";
      }
      else 
      {
        std::vector<Vec3> newStreamline;
        int presumedLength = int(currentStreamline.size()*targetElementLength_/this->lineStepWidth_+10);
        newStreamline.reserve(presumedLength);
        
        VLOG(1) << "streamline no " << i << ", reserve length " << presumedLength;
        VLOG(1) << "targetElementLength_: " << targetElementLength_ << ", scalingFactor: " << scalingFactor;
        
        Vec3 lastPoint = currentStreamline.front()*scalingFactor;
        Vec3 previousStreamlinePoint = lastPoint;  // last point that was inserted into the new streamline
        // use starting point of streamline
        newStreamline.push_back(lastPoint);
        double length = 0.0;
        bool firstPoint = true;
        
        // loop over points of streamline
        for (std::vector<Vec3>::iterator pointsIter = currentStreamline.begin(); pointsIter != currentStreamline.end(); pointsIter++)
        {
          if (!firstPoint)
          {
            Vec3 currentPoint = (*pointsIter)*scalingFactor;
            // sum up length since last element started
            length += MathUtility::length<3>(currentPoint - lastPoint);
            
            VLOG(1) << "old streamline interval " << lastPoint << " - " << currentPoint << ", new lentgh: " << length << " (targetElementLength=" << targetElementLength_ << ")";
            
            if (length > targetElementLength_)
            {
              double alpha = targetElementLength_/length;
              Vec3 point = (1. - alpha) * previousStreamlinePoint + alpha * currentPoint;
             
              VLOG(1) << "  length is too big, alpha=" << alpha << ", take intermediate point " << point
                << ", distance to previous point " << previousStreamlinePoint << ": " << MathUtility::length<3>(point - previousStreamlinePoint);
              
              newStreamline.push_back(point);
              
              previousStreamlinePoint = point;
              lastPoint = previousStreamlinePoint;
              length = 0.0;
            }
            else 
            {
              lastPoint = currentPoint;
            }
            
          }
          firstPoint = false;
        }
        LOG(DEBUG) << "Scaled streamline by factor " << scalingFactor << ", resampled from lineStepWidth " << this->lineStepWidth_ << " to targetElementLength " << targetElementLength_
          << ", now it has " << newStreamline.size() << " points, length " << lengths[i]*scalingFactor;
        streamlines[i] = newStreamline;
      }
    }
  }
    
  LOG(DEBUG) << "Number of streamlines after resampling: " << streamlines.size();
}


};
