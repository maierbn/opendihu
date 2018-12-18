#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
traceStreamlines(int nRanksZ, int rankZNo, double streamlineDirection, bool streamlineDirectionUpwards, std::vector<Vec3> &seedPoints, std::vector<std::vector<Vec3>> &streamlinePoints, const std::array<std::shared_ptr<FunctionSpaceType>,4> &ghostMesh)
{
  int nStreamlines = seedPoints.size();
  streamlinePoints.resize(nStreamlines);

  if (nRanksZ == 1)
  {
    // if there is only one rank, trace streamline from the center in both directions
    for (int face = (int)Mesh::face_t::face0Minus; face <= (int)Mesh::face_t::face1Plus; face++)
    {
      this->functionSpace_->setGhostMesh(Mesh::face_t::face0Minus, ghostMesh[face]);
    }

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
      std::stringstream name;
      name << "04_streamline_" << i << "_";
      PyObject_CallFunction(functionOutputPoints_, "s i O f", name.str().c_str(), currentRankSubset_->ownRankNo(),
                            PythonUtility::convertToPython<std::vector<Vec3>>::get(streamlinePoints[i]), 0.5);
      PythonUtility::checkForError();
#endif
    }
  }
  else
  {
    // multiple ranks

    // determine if previously set seedPoints are used or if they are received from neighbouring rank
    if (nRanksZ > 1 && rankZNo != int(nRanksZ/2) && rankZNo != int(nRanksZ/2+1))
    {
      int neighbourRankNo;
      if (streamlineDirectionUpwards)
      {
        neighbourRankNo = meshPartition_->neighbourRank(Mesh::face_t::face2Minus);
      }
      else
      {
        neighbourRankNo = meshPartition_->neighbourRank(Mesh::face_t::face2Plus);
      }

      // receive seed points
      MPIUtility::handleReturnValue(MPI_Recv(seedPoints.data(), seedPoints.size(), MPI_DOUBLE, neighbourRankNo,
                                            0, currentRankSubset_->mpiCommunicator(), MPI_STATUS_IGNORE), "MPI_Recv");
    }

    LOG(DEBUG) << " on " << nRanksZ << " ranks in Z direction, trace " << nStreamlines << " streamlines";

    // trace streamlines from seed points
    for (int i = 0; i < nStreamlines; i++)
    {
      Vec3 &startingPoint = seedPoints[i];
      streamlinePoints[i].push_back(startingPoint);

      for (int face = (int)Mesh::face_t::face0Minus; face <= (int)Mesh::face_t::face1Plus; face++)
      {
        this->functionSpace_->setGhostMesh(Mesh::face_t::face0Minus, ghostMesh[face]);
      }

      this->traceStreamline(startingPoint, streamlineDirection, streamlinePoints[i]);

#ifndef NDEBUG
      std::stringstream name;
      name << "04_streamline_" << i << "_";
      PyObject_CallFunction(functionOutputPoints_, "s i O f", name.str().c_str(), currentRankSubset_->ownRankNo(),
                            PythonUtility::convertToPython<std::vector<Vec3>>::get(streamlinePoints[i]), 1.0);
      PythonUtility::checkForError();
#endif

    }

    // send end points of streamlines to next rank that continues the streamline
    if (nRanksZ > 1 && rankZNo != nRanksZ-1 && rankZNo != 0)
    {
      // fill send buffer
      std::vector<double> sendBuffer(nStreamlines);
      for (int streamlineIndex = 0; streamlineIndex < nStreamlines; streamlineIndex++)
      {
        for (int i = 0; i < 3; i++)
        {
          sendBuffer[streamlineIndex] = streamlinePoints[streamlineIndex].back()[i];
        }
      }

      int neighbourRankNo;
      if (streamlineDirectionUpwards)
      {
        neighbourRankNo = meshPartition_->neighbourRank(Mesh::face_t::face2Plus);
      }
      else
      {
        neighbourRankNo = meshPartition_->neighbourRank(Mesh::face_t::face2Minus);
      }

      // send end points of streamlines
      MPIUtility::handleReturnValue(MPI_Send(sendBuffer.data(), sendBuffer.size(), MPI_DOUBLE, neighbourRankNo,
                                            0, currentRankSubset_->mpiCommunicator()), "MPI_Send");
    }

    // reorder streamline points such that they go from bottom to top
    if (streamlineDirection < 0)
    {
      for (int streamlineIndex = 0; streamlineIndex < nStreamlines; streamlineIndex++)
      {
        std::reverse(streamlinePoints[streamlineIndex].begin(), streamlinePoints[streamlineIndex].end());
      }
    }
  }
}

};  // namespace
