#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
traceStreamlines(int nRanksZ, int rankZNo, double streamlineDirection, bool streamlineDirectionUpwards, std::vector<Vec3> &seedPoints, std::vector<std::vector<Vec3>> &streamlinePoints)
{
  int nStreamlines = seedPoints.size();
  streamlinePoints.resize(nStreamlines);

  if (nRanksZ == 1)
  {
    // trace streamlines from seed points
    int nStreamlines = seedPoints.size();
    LOG(DEBUG) << " on 1 rank, trace " << nStreamlines << " streamlines";
    for (int i = 0; i < nStreamlines; i++)
    {
      // get starting point
      Vec3 &startingPoint = seedPoints[i];

      // trace streamlines forwards
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
      streamlinePoints[i].insert(streamlinePoints[i].begin(), backwardPoints.rbegin(), backwardPoints.rend());
      streamlinePoints[i].insert(streamlinePoints[i].end(), startingPoint);
      streamlinePoints[i].insert(streamlinePoints[i].end(), forwardPoints.begin(), forwardPoints.end());

      // now streamline is in order from bottom (Z-) to top (Z+)
      LOG(DEBUG) << "i=" << i << ", n points in streamline: " << streamlinePoints[i].size();
      VLOG(1) << "streamline: " << streamlinePoints[i];

#ifndef NDEBUG
#ifdef STL_OUTPUT
#ifdef STL_OUTPUT_VERBOSE
      std::stringstream name;
      name << "04_streamline_" << i << "_";
      PyObject_CallFunction(functionOutputPoints_, "s i i O f", name.str().c_str(), currentRankSubset_->ownRankNo(), level_,
                            PythonUtility::convertToPython<std::vector<Vec3>>::get(streamlinePoints[i]), 0.5);
      PythonUtility::checkForError();
#endif
#endif
#endif
    }
  }
  else
  {
    // multiple ranks

    // The overall picture is that global streamlines begin at the center (at rankZNo/2).
    // Rank int(nRanksZ/2) send the initial seed points to the rank below (int(nRanksZ/2)-1)
    // Then every rank traces its streamlines and sends the end points as new seed points to the next rank (lower or upper neighbour, depending on streamlineDirection)

    LOG(DEBUG) << "call exchangeBorderSeedPoints with " << seedPoints.size() << " seed points";

    // determine if previously set seedPoints are used or if they are received from neighbouring rank, receive seed points or send them to lower neighbour, if own rank is int(nRanksZ/2)
    exchangeBorderSeedPoints(nRanksZ, rankZNo, streamlineDirectionUpwards, seedPoints);

    // receive seed points from neighbouring rank
    exchangeBorderSeedPointsBeforeTracing(nRanksZ, rankZNo, streamlineDirectionUpwards, seedPoints);

    LOG(DEBUG) << " on " << nRanksZ << " ranks in Z direction, trace " << nStreamlines << " streamlines, streamlineDirectionUpwards: " << streamlineDirectionUpwards;

#ifndef NDEBUG
#ifdef STL_OUTPUT
    PyObject_CallFunction(functionOutputPoints_, "s i i O f", "03_seed_points", currentRankSubset_->ownRankNo(), level_,
                          PythonUtility::convertToPython<std::vector<Vec3>>::get(seedPoints), 0.2);
    PythonUtility::checkForError();
#endif
#endif

    //MPI_Barrier(currentRankSubset_->mpiCommunicator());
    //LOG(FATAL) << "end before tracing of streamlines";

    // trace streamlines from seed points
    for (int i = 0; i < nStreamlines; i++)
    {
      Vec3 &startingPoint = seedPoints[i];
      streamlinePoints[i].push_back(startingPoint);


      /* debugging condition, TODO: remove */
      /*int ownRankNo = currentRankSubset_->ownRankNo();
      if (!((ownRankNo == 1)))
      {
        continue;
      }*/
      /* end */

      this->traceStreamline(startingPoint, streamlineDirection, streamlinePoints[i]);

      // if everything was cleared, add seed point
      if (streamlinePoints[i].empty())
      {
        streamlinePoints[i].push_back(startingPoint);


#ifndef NDEBUG
#ifdef STL_OUTPUT
//#ifdef STL_OUTPUT_VERBOSE
        std::stringstream name;
        name << "04_raw_invalid_streamline_" << i << "_";
        PyObject_CallFunction(functionOutputStreamline_, "s i i O f", name.str().c_str(), currentRankSubset_->ownRankNo(), level_,
                              PythonUtility::convertToPython<std::vector<Vec3>>::get(streamlinePoints[i]), 0.1);
        PythonUtility::checkForError();
//#endif
#endif
#endif
      }
    }

#ifndef NDEBUG
#ifdef STL_OUTPUT
      std::stringstream name;
      name << "04_raw_streamlines_";
      PyObject_CallFunction(functionOutputStreamlines_, "s i i O f", name.str().c_str(), currentRankSubset_->ownRankNo(), level_,
                            PythonUtility::convertToPython<std::vector<std::vector<Vec3>>>::get(streamlinePoints), 0.1);
      PythonUtility::checkForError();
#endif
#endif

    // reorder streamline points such that they go from bottom to top
    if (streamlineDirection < 0)
    {
      VLOG(1) << "reverse streamlines direction";
      for (int streamlineIndex = 0; streamlineIndex < nStreamlines; streamlineIndex++)
      {
        //VLOG(1) << streamlineIndex << " has " << streamlinePoints[streamlineIndex].size() << " points, before: " << streamlinePoints[streamlineIndex];
        std::reverse(streamlinePoints[streamlineIndex].begin(), streamlinePoints[streamlineIndex].end());
        //LOG(DEBUG) << streamlineIndex << " after: " << streamlinePoints[streamlineIndex];
      }
    }

    LOG(DEBUG) << "call exchangeSeedPointsAfterTracing with " << seedPoints.size() << " seed points from traceStreamlines";

    // send end points of streamlines to next rank that continues the streamline
    exchangeBorderSeedPointsAfterTracing(nRanksZ, rankZNo, streamlineDirectionUpwards, streamlinePoints);

  }

  //MPI_Barrier(currentRankSubset_->mpiCommunicator());
  //LOG(FATAL) << "end after all streamlines were traced";
}

} // namespace
