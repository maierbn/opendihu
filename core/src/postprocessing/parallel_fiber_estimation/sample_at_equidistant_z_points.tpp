#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
sampleAtEquidistantZPoints(std::vector<std::vector<Vec3>> &streamlinePoints, const std::vector<Vec3> &seedPoints, std::vector<std::vector<Vec3>> &streamlineZPoints)
{
  // determine z range of current subdomain
  int nRanksZ = meshPartition_->nRanks(2);
  int rankZNo = meshPartition_->ownRankPartitioningIndex(2);

  double zRangeTotal = topZClip_ - bottomZClip_;
  double zRangeCurrentLevel = zRangeTotal / nRanksZ;
  double bottomZClip = bottomZClip_ + zRangeCurrentLevel*rankZNo;
  double topZClip = bottomZClip_ + zRangeCurrentLevel*(rankZNo+1);

  LOG(DEBUG) << "z bounds: " << bottomZClip << ", " << topZClip;

  int nStreamlines = streamlinePoints.size();
  streamlineZPoints.resize(nStreamlines);

  // loop over all traced streamlines in this subdomain
  for (int i = 0; i < nStreamlines; i++)
  {
    LOG(DEBUG) << " streamline " << i << " has " << streamlinePoints[i].size() << " points.";

    // the streamline is expected to have at least one point, the seed point
    assert(!streamlinePoints[i].empty());

    if (streamlinePoints[i].size() == 1)
    {
      streamlineZPoints[i].push_back(streamlinePoints[i][0]);
      continue;
    }

    // here streamlinePoints contains at least 2 points
    std::vector<Vec3>::const_iterator streamlineIter = streamlinePoints[i].begin();

    // loop over z levels
    double currentZ;
    for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZNew_; zLevelIndex++)
    {
      currentZ = bottomZClip + double(zLevelIndex) / (nBorderPointsZNew_-1) * (topZClip - bottomZClip);
      VLOG(1) << "currentZ: " << currentZ;

      while(streamlineIter != streamlinePoints[i].end())
      {
        VLOG(1) << "  z: " << (*streamlineIter)[2];
        if ((*streamlineIter)[2] > currentZ + 1e-9)
          break;

        streamlineIter++;
        if (streamlineIter == streamlinePoints[i].end())
          break;
      }

      if (streamlineIter == streamlinePoints[i].end())
      {
        VLOG(1) << "end";
        streamlineIter--;
      }

      Vec3 currentPoint = *streamlineIter;
      Vec3 previousPoint = currentPoint;

      // if there was no point before the current z, this means the streamline begins way higher than the current z
      if (streamlineIter == streamlinePoints[i].begin())
      {
        VLOG(1) << "no previous point, skip. first currentPoint: " << currentPoint;
        continue;
      }
      else
      {
        previousPoint = *(streamlineIter-1);
      }

      VLOG(1) << "currentPoint: " << currentPoint << ", previousPoint: " << previousPoint;

      // now previousPoint is the last point under currentZ and currentPoint is the first over currentZ
      if (fabs(currentPoint[2] - previousPoint[2]) < 1e-9)
      {
        VLOG(1) << "same, continue";
        continue;
      }

      double alpha = (currentZ - previousPoint[2]) / (currentPoint[2] - previousPoint[2]);
      Vec3 point = (1.-alpha) * previousPoint + alpha * currentPoint;
      streamlineZPoints[i].push_back(point);
      VLOG(1) << "alpha: " << alpha;
    }

    LOG(DEBUG) << " n sampled points: " << streamlineZPoints[i].size() << ", nBorderPointsXNew_: " << nBorderPointsXNew_ << ", nBorderPointsZNew_: " << nBorderPointsZNew_;

#ifndef NDEBUG
#ifdef STL_OUTPUT
//#ifdef STL_OUTPUT_VERBOSE
    std::stringstream name;
    name << "05_sampled_streamline_" << i << "_";
    PyObject_CallFunction(functionOutputStreamline_, "s i O f", name.str().c_str(), currentRankSubset_->ownRankNo(),
                          PythonUtility::convertToPython<std::vector<Vec3>>::get(streamlineZPoints[i]), 0.1);
    PythonUtility::checkForError();
//#endif
#endif
#endif

    // if streamline is not complete
    if (streamlineZPoints[i].size() != nBorderPointsZNew_)
    {
      LOG(ERROR) << "Streamline " << i << " is not complete, i.e. does not run from \"bottomZClip\" to \"topZClip\" .";

      // assign seed point instead of incomplete streamline
      streamlineZPoints[i].resize(1);
      streamlineZPoints[i][0] = seedPoints[i];
    }
  }

}

};  // namespace
