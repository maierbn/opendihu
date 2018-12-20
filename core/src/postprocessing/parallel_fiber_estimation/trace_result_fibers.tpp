#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
bool ParallelFiberEstimation<BasisFunctionType>::
checkTraceFinalFibers(int &level)
{
  // determine current level = log2(nRanksPerCoordinateDirection_)
  level = 0;
  int nRanksPerCoordinateDirection = nRanksPerCoordinateDirection_[0];
  while (nRanksPerCoordinateDirection >>= 1)
  {
    level++;
  }

  int nRanksAvailable = this->context_.partitionManager()->nRanksCommWorld();

  // decide if the algorithm should no more refine the subdomains but trace the final fibers
  bool traceFinalFibers = false;
  if (level == maxLevel_)
  {
    traceFinalFibers = true;
  }

  if (!traceFinalFibers)
  {
    if (nRanksAvailable < currentRankSubset_->size()*8)
    {
      LOG(WARNING) << "Cannot run algorithm until level " << maxLevel_ << ", currently at level " << level
        << ", total number of ranks is " << nRanksAvailable << ", number needed for next level would be " << currentRankSubset_->size()*8 << ".";
      //traceFinalFibers = true;
    }
  }
  return traceFinalFibers;
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
traceResultFibers(double streamlineDirection, int seedPointsZIndex, const std::vector<Vec3> &nodePositions)
{
  LOG(DEBUG) << "final level";

  int subdomainNNodesX = nBorderPointsXNew_;
  int subdomainNNodesY = nBorderPointsXNew_;
  std::vector<Vec3> seedPoints;

  // create seed points
  for (int j = 0; j < nBorderPointsXNew_-1; j++)
  {
    for (int i = 0; i < nBorderPointsXNew_-1; i++)
    {
      Vec3 p0 = nodePositions[seedPointsZIndex*subdomainNNodesX*subdomainNNodesY + j*subdomainNNodesX + i];
      Vec3 p1 = nodePositions[seedPointsZIndex*subdomainNNodesX*subdomainNNodesY + j*subdomainNNodesX + i+1];
      Vec3 p2 = nodePositions[seedPointsZIndex*subdomainNNodesX*subdomainNNodesY + (j+1)*subdomainNNodesX + i];
      Vec3 p3 = nodePositions[seedPointsZIndex*subdomainNNodesX*subdomainNNodesY + (j+1)*subdomainNNodesX + i+1];
      seedPoints.push_back(0.25*(p0 + p1 + p2 + p3));
    }
  }

#ifndef NDEBUG
  PyObject_CallFunction(functionOutputPoints_, "s i O f", "03_final_seed_points", currentRankSubset_->ownRankNo(),
                        PythonUtility::convertToPython<std::vector<Vec3>>::get(seedPoints), 1.0);
  PythonUtility::checkForError();
#endif

  // trace streamlines from seed points
  int nStreamlines = seedPoints.size();
  std::vector<std::vector<Vec3>> streamlinePoints(nStreamlines);
  for (int i = 0; i < nStreamlines; i++)
  {
    Vec3 &startingPoint = seedPoints[i];
    streamlinePoints[i].push_back(startingPoint);

    this->traceStreamline(startingPoint, streamlineDirection, streamlinePoints[i]);
  }

  // exchange end points of streamlines and adjust them
}

};  // namespace
