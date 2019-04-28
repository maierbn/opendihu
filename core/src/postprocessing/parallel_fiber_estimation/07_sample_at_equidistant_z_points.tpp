#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
sampleAtEquidistantZPoints(std::vector<std::vector<Vec3>> &streamlinePoints, const std::vector<Vec3> &seedPoints, std::vector<std::vector<Vec3>> &streamlineZPoints)
{
  // determine z range of current subdomain
  double bottomZClip = 0;
  double topZClip = 0;
  computeBottomTopZClip(bottomZClip, topZClip);

  LOG(DEBUG) << "z bounds: " << bottomZClip << ", " << topZClip;

  int nStreamlines = streamlinePoints.size();
  streamlineZPoints.resize(nStreamlines);

  // loop over all traced streamlines in this subdomain
  for (int i = 0; i < nStreamlines; i++)
  {
    LOG(DEBUG) << " streamline " << i << " has " << streamlinePoints[i].size() << " points.";

    sampleStreamlineAtEquidistantZPoints(streamlinePoints[i], seedPoints[i], bottomZClip, topZClip, streamlineZPoints[i], i);
  }
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
computeBottomTopZClip(double &bottomZClip, double &topZClip)
{
  // determine z range of current subdomain
  int nRanksZ = meshPartition_->nRanks(2);
  int rankZNo = meshPartition_->ownRankPartitioningIndex(2);

  double zRangeTotal = topZClip_ - bottomZClip_;
  double zRangeCurrentLevel = zRangeTotal / nRanksZ;
  bottomZClip = bottomZClip_ + zRangeCurrentLevel*rankZNo;
  topZClip = bottomZClip_ + zRangeCurrentLevel*(rankZNo+1);
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
sampleStreamlineAtEquidistantZPoints(std::vector<Vec3> &streamlinePoints, const Vec3 &seedPoint, double bottomZClip, double topZClip,
                                     std::vector<Vec3> &streamlineZPoints, int streamlineNoDebugging)
{
  // the streamline is expected to have at least one point, the seed point
  assert(!streamlinePoints.empty());
  streamlineZPoints.clear();
  LOG(DEBUG) << "sampleStreamlineAtEquidistantZPoints, streamlinePoints: " << streamlinePoints.size()
    << ", first: " << streamlinePoints[0] << ", last: " << streamlinePoints[streamlinePoints.size()-1];

  if (streamlinePoints.size() == 1)
  {
    streamlineZPoints.push_back(streamlinePoints[0]);
    return;
  }

  // here streamlinePoints contains at least 2 points
  std::vector<Vec3>::const_iterator streamlineIter = streamlinePoints.begin();

  // loop over z levels of the streamline
  const double epsilon = 1e-5;   // this has to be ~1e-5 because the seed points are perturbed by MPI communcation around this amount
  double currentZ;
  for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZNew_; zLevelIndex++)
  {
    // compute current z level at which a point is searched
    currentZ = bottomZClip + double(zLevelIndex) / (nBorderPointsZNew_-1) * (topZClip - bottomZClip);
    VLOG(1) << "zLevelIndex " << zLevelIndex << ", currentZ: " << currentZ;

    // advance streamline until current z is reached
    while(streamlineIter != streamlinePoints.end())
    {
      VLOG(1) << "  z: " << (*streamlineIter)[2]  << " (float: " << (*streamlineIter)[2] - int((*streamlineIter)[2])
        << "), currentZ+epsilon = " << currentZ + epsilon << " (float: " << (currentZ + epsilon) - int(currentZ + epsilon) << ") "
        << "diff: " << ((*streamlineIter)[2] - (currentZ + epsilon)) << ", higher: " << std::boolalpha << ((*streamlineIter)[2] > currentZ + epsilon);
      if ((*streamlineIter)[2] > currentZ + epsilon)
        break;

      VLOG(1) << "  streamlineIter++";
      streamlineIter++;
      if (streamlineIter == streamlinePoints.end())
        break;
    }

    if (streamlineIter == streamlinePoints.end())
    {
      VLOG(1) << "end";
      streamlineIter--;

      if (topZClip - (*streamlineIter)[2] > (topZClip - bottomZClip)*0.05)
      {
        LOG(DEBUG) << "Streamline does not reach topZClip: " << topZClip << ", but finishes at " << (*streamlineIter)[2]
          << " (" << (topZClip - (*streamlineIter)[2]) / (topZClip - bottomZClip) * 100.0 << "% too early)";
        break;
      }
    }

    Vec3 currentPoint = *streamlineIter;
    Vec3 previousPoint = currentPoint;

    // if there was no point before the current z, this means the streamline begins way higher than the current z
    if (streamlineIter == streamlinePoints.begin())
    {
      if ((*streamlineIter)[2] < currentZ + 0.3 && streamlineIter+1 != streamlinePoints.end())
      {
        // if the streamline begins a little (0.3) higher than the bottom clip but not too much
        currentPoint = *(streamlineIter+1);
      }
      else
      {
        // if the streamline begins way higher
        VLOG(1) << "no previous point, skip. first currentPoint: " << currentPoint;
        continue;
      }
    }
    else
    {
      VLOG(1) << "previous point is streamlineIter-1";
      previousPoint = *(streamlineIter-1);
    }

    VLOG(1) << "currentPoint: " << currentPoint << ", previousPoint: " << previousPoint;

    // now previousPoint is the last point under currentZ and currentPoint is the first over currentZ

    // if point is a duplicate, skip
    if (fabs(currentPoint[2] - previousPoint[2]) < epsilon)
    {
      VLOG(1) << "same, continue";
      continue;
    }

    double alpha = (currentZ - previousPoint[2]) / (currentPoint[2] - previousPoint[2]);
    Vec3 point = (1.-alpha) * previousPoint + alpha * currentPoint;
    streamlineZPoints.push_back(point);
    VLOG(1) << "alpha: " << alpha << ", take point " << point << ", streamline now has " << streamlineZPoints.size() << " points";
  }

  LOG(DEBUG) << " n sampled points: " << streamlineZPoints.size()
    << ", clip: [" << bottomZClip << "," << topZClip << "], first: " << streamlineZPoints[0]
    << ", last: " << streamlineZPoints[streamlineZPoints.size()-1]
    << ", nBorderPointsXNew_: " << nBorderPointsXNew_ << ", nBorderPointsZNew_: " << nBorderPointsZNew_;

#ifndef NDEBUG
#ifdef STL_OUTPUT
//#ifdef STL_OUTPUT_VERBOSE
  std::stringstream name;
  name << "05_sampled_streamline_" << streamlineNoDebugging << "_";
  PyObject_CallFunction(functionOutputStreamline_, "s i i O f", name.str().c_str(), currentRankSubset_->ownRankNo(), level_,
                        PythonUtility::convertToPython<std::vector<Vec3>>::get(streamlineZPoints), 0.1);
  PythonUtility::checkForError();
//#endif
#endif
#endif

  // if streamline is not complete
  if (streamlineZPoints.size() != nBorderPointsZNew_)
  {
    LOG(DEBUG) << "Streamline " << streamlineNoDebugging << " is not complete, i.e. does not run from \"bottomZClip\" to \"topZClip\" .";

    // assign seed point instead of incomplete streamline
    streamlineZPoints.resize(1);
    streamlineZPoints[0] = streamlinePoints[0];
  }

}

} // namespace
