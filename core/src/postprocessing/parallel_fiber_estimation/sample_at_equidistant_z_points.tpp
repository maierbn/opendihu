#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
sampleAtEquidistantZPoints(std::vector<std::vector<Vec3>> &streamlinePoints, std::vector<std::vector<Vec3>> &streamlineZPoints)
{
  double currentZ = bottomZClip_;
  double zIncrement = double(topZClip_ - bottomZClip_) / (nBorderPointsZNew_ - 1.0);

  LOG(DEBUG) << "z bounds: " << bottomZClip_ << ", " << topZClip_ << ", zIncrement: " << zIncrement;

  int nStreamlines = streamlinePoints.size();
  streamlineZPoints.resize(nStreamlines);
  for (int i = 0; i < nStreamlines; i++)
  {
    LOG(DEBUG) << " streamline " << i << " has " << streamlinePoints[i].size() << " points.";

    // the streamline is expected to have at least one point, the seed point
    assert(!streamlinePoints[i].empty());

    Vec3 previousPoint = streamlinePoints[i].front();
    streamlineZPoints[i].reserve(nBorderPointsZNew_);
    streamlineZPoints[i].push_back(previousPoint);

    if (streamlinePoints[i].size() > 1)
    {
      assert(streamlinePoints[i].size() > 1);

      currentZ = bottomZClip_ + zIncrement;
      LOG(DEBUG) << "first point: " << streamlinePoints[i].front() << ", currentZ: " << currentZ;

      for (std::vector<Vec3>::const_iterator iter = streamlinePoints[i].begin()+1; iter != streamlinePoints[i].end(); iter++)
      {
        const Vec3 &currentPoint = *iter;
        if (currentPoint[2] > currentZ || iter+1 == streamlinePoints[i].end())
        {
          double alpha = (currentZ - previousPoint[2]) / (currentPoint[2] - previousPoint[2]);
          Vec3 point = (1.-alpha) * previousPoint + alpha * currentPoint;
          streamlineZPoints[i].push_back(point);

          currentZ += zIncrement;
        }
        previousPoint = currentPoint;
      }

  #ifndef NDEBUG
      std::stringstream name;
      name << "05_sampled_streamline_" << i << "_";
      PyObject_CallFunction(functionOutputPoints_, "s i O f", name.str().c_str(), currentRankSubset_->ownRankNo(),
                            PythonUtility::convertToPython<std::vector<Vec3>>::get(streamlineZPoints[i]), 1.0);
      PythonUtility::checkForError();
  #endif

      LOG(DEBUG) << "n sampled points: " << streamlineZPoints[i].size() << ", nBorderPointsXNew_: " << nBorderPointsXNew_ << ", nBorderPointsZNew_: " << nBorderPointsZNew_;
    }

    if (streamlineZPoints[i].size() != nBorderPointsZNew_)
    {
      LOG(ERROR) << "Streamline " << i << " is not complete, i.e. does not run from \"bottomZClip\" to \"topZClip\" "
        << std::endl << "Try adjusting \"bottomZClip\" and \"topZClip\" or mesh width.";
      streamlineZPoints[i].resize(1);
    }
    //assert(streamlineZPoints[i].size() == nBorderPointsZNew_);
  }
}

};  // namespace
