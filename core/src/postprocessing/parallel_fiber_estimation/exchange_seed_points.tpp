#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
exchangeSeedPointsBeforeTracing(int nRanksZ, int rankZNo, bool streamlineDirectionUpwards, std::vector<Vec3> &seedPoints)
{
  // determine if previously set seedPoints are used or if they are received from neighbouring rank
  LOG(DEBUG) << "exchangeSeedPointsBeforeTracing, rankZNo: " << rankZNo << ", streamlineDirectionUpwards: " << streamlineDirectionUpwards << ", int(nRanksZ/2): " << int(nRanksZ/2);

  if (nRanksZ == 1)
    return;

  if (rankZNo != int(nRanksZ/2))
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
    int tag = currentRankSubset_->ownRankNo()*100 + neighbourRankNo*10000 + level_*10 + 7;
    if (rankZNo == int(nRanksZ/2)-1)
    {
      tag = currentRankSubset_->ownRankNo()*100 + neighbourRankNo*10000 + level_*10 + 6;
    }

    LOG(DEBUG) << "receive " << seedPoints.size()*3 << " seed point values (" << seedPoints.size() << " seed points) from rank " << neighbourRankNo << " (tag: " << tag << ")";
    std::vector<double> receiveBuffer(seedPoints.size()*3);
    MPIUtility::handleReturnValue(MPI_Recv(receiveBuffer.data(), receiveBuffer.size(), MPI_DOUBLE, neighbourRankNo,
                                           tag, currentRankSubset_->mpiCommunicator(), MPI_STATUS_IGNORE), "MPI_Recv");

    // fill seed points from receive buffer
    for (int seedPointIndex = 0; seedPointIndex < seedPoints.size(); seedPointIndex++)
    {
      for (int i = 0; i < 3; i++)
      {
        seedPoints[seedPointIndex][i] = receiveBuffer[seedPointIndex*3 + i];
      }
    }
    //LOG(DEBUG) << "received " << seedPoints.size() << " seed points from rank " << neighbourRankNo << ": " << seedPoints;

  }

  // on rank int(nRanksZ/2), send seed points to rank below
  //  rank int(nRanksZ/2)     |_|_| tracing direction: ^
  //  rank int(nRanksZ/2)-1   | | | tracing direction: v
  if (rankZNo == int(nRanksZ/2))
  {
    int neighbourRankNo = meshPartition_->neighbourRank(Mesh::face_t::face2Minus);

    // fill send buffer
    std::vector<double> sendBuffer(seedPoints.size()*3);
    for (int seedPointIndex = 0; seedPointIndex < seedPoints.size(); seedPointIndex++)
    {
      for (int i = 0; i < 3; i++)
      {
        sendBuffer[seedPointIndex*3 + i] = seedPoints[seedPointIndex][i];
      }
    }


    // send seed points
    int tag = currentRankSubset_->ownRankNo()*10000 + neighbourRankNo*100 + level_*10 + 6;
    LOG(DEBUG) << "send " << sendBuffer.size() << " initial seed points values to rank " << neighbourRankNo << " (tag: " << tag << ")";
    MPIUtility::handleReturnValue(MPI_Send(sendBuffer.data(), sendBuffer.size(), MPI_DOUBLE, neighbourRankNo,
                                            tag, currentRankSubset_->mpiCommunicator()), "MPI_Send");
  }
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
exchangeSeedPointsAfterTracing(int nRanksZ, int rankZNo, bool streamlineDirectionUpwards, std::vector<Vec3> &seedPoints, std::vector<std::vector<Vec3>> &streamlinePoints)
{
  LOG(DEBUG) << "exchangeSeedPointsAfterTracing, nRanksZ: " << nRanksZ << ", rankZNo: " << rankZNo << ", streamlineDirectionUpwards: " << streamlineDirectionUpwards
    << ", ownRankPartitioningIndex_: " << this->meshPartition_->ownRankPartitioningIndex(0) << "," << this->meshPartition_->ownRankPartitioningIndex(1) << "," << this->meshPartition_->ownRankPartitioningIndex(2);

  if (nRanksZ == 1)
    return;

  // this assumes the streamline always goes from bottom to top

  // send end points of streamlines to next rank that continues the streamline
  if (nRanksZ > 1 && rankZNo != nRanksZ-1 && rankZNo != 0)
  {

    int neighbourRankNo;
    if (streamlineDirectionUpwards)
    {
      neighbourRankNo = meshPartition_->neighbourRank(Mesh::face_t::face2Plus);
    }
    else
    {
      neighbourRankNo = meshPartition_->neighbourRank(Mesh::face_t::face2Minus);
    }

    // fill send buffer
    std::vector<double> sendBuffer(seedPoints.size()*3);
    for (int streamlineIndex = 0; streamlineIndex < seedPoints.size(); streamlineIndex++)
    {
      int streamlinePointNo = 0;

      // if the streamline is going upwards, the next seed point is the upper most, i.e. the last, otherwise it is the first
      if (streamlineDirectionUpwards)
      {
        streamlinePointNo = streamlinePoints[streamlineIndex].size();
      }

      for (int i = 0; i < 3; i++)
      {
        sendBuffer[streamlineIndex*3+i] = streamlinePoints[streamlineIndex][streamlinePointNo][i];
      }
    }

    int tag = currentRankSubset_->ownRankNo()*10000 + neighbourRankNo*100 + level_*10 + 7;
    LOG(DEBUG) << "send " << sendBuffer.size() << " seed point values to " << neighbourRankNo << ", tag: " << tag;
    // send end points of streamlines
    MPIUtility::handleReturnValue(MPI_Send(sendBuffer.data(), sendBuffer.size(), MPI_DOUBLE, neighbourRankNo,
                                          tag, currentRankSubset_->mpiCommunicator()), "MPI_Send");
  }
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
exchangeSeedPointsAfterTracingKeyFibers(int nRanksZ, int rankZNo, bool streamlineDirectionUpwards, int nFibersX, std::vector<Vec3> &seedPoints, std::vector<std::vector<Vec3>> &streamlinePoints)
{
  // nRanksZ, rankZNo, streamlineDirectionUpwards, nFibersX, seedPoints, fibers
  // seedPoints size is nBorderPointsXNew_ * nBorderPointsXNew_
  // streamlinePoints size is nFibersX * nFibersX

  LOG(DEBUG) << "exchangeSeedPointsAfterTracingKeyFibers, nRanksZ: " << nRanksZ << ", rankZNo: " << rankZNo << ", streamlineDirectionUpwards: " << streamlineDirectionUpwards
    << ", ownRankPartitioningIndex_: " << this->meshPartition_->ownRankPartitioningIndex(0) << "," << this->meshPartition_->ownRankPartitioningIndex(1) << "," << this->meshPartition_->ownRankPartitioningIndex(2);
  LOG(DEBUG) << "n seed points: " << seedPoints.size() << ", n streamlines: " << streamlinePoints.size();

  std::stringstream stream;
  for (int i = 0; i < streamlinePoints.size(); i++)
  {
    stream << " " << streamlinePoints[i].size();
  }
  LOG(DEBUG) << " streamline sizes: " << stream.str();

  if (nRanksZ == 1)
    return;

  // this assumes the streamline always goes from bottom to top

  // send end points of streamlines to next rank that continues the streamline
  if (nRanksZ > 1 && rankZNo != nRanksZ-1 && rankZNo != 0)
  {

    int neighbourRankNo;
    if (streamlineDirectionUpwards)
    {
      neighbourRankNo = meshPartition_->neighbourRank(Mesh::face_t::face2Plus);
      LOG(DEBUG) << "get 2+ neighbourRankNo: " << neighbourRankNo;
    }
    else
    {
      neighbourRankNo = meshPartition_->neighbourRank(Mesh::face_t::face2Minus);
      LOG(DEBUG) << "get 2- neighbourRankNo: " << neighbourRankNo;
    }

    // fill send buffer
    std::vector<double> sendBuffer(seedPoints.size()*3);
      
    for (int j = 0; j < nBorderPointsXNew_; j++)
    {
      for (int i = 0; i < nBorderPointsXNew_; i++)
      {
        int streamlineIndex = j * (nFineGridFibers_+1) * nFibersX + i * (nFineGridFibers_+1);

        int streamlinePointNo = 0;

        // if the streamline is going upwards, the next seed point is the upper most, i.e. the last, otherwise it is the first
        if (streamlineDirectionUpwards)
        {
          streamlinePointNo = streamlinePoints[streamlineIndex].size() - 1;
        }

        for (int k = 0; k < 3; k++)
        {
          sendBuffer[streamlineIndex*3+k] = streamlinePoints[streamlineIndex][streamlinePointNo][k];
        }
      }
    }

    int tag = currentRankSubset_->ownRankNo()*10000 + neighbourRankNo*100 + level_*10 + 7;
    LOG(DEBUG) << "send " << sendBuffer.size() << " seed point values (" << sendBuffer.size()/3 << " seed points) to " << neighbourRankNo << ", tag: " << tag;
    LOG(DEBUG) << "sendBuffer: " << sendBuffer;

    // send end points of streamlines
    MPIUtility::handleReturnValue(MPI_Send(sendBuffer.data(), sendBuffer.size(), MPI_DOUBLE, neighbourRankNo,
                                          tag, currentRankSubset_->mpiCommunicator()), "MPI_Send");
  }
}


} // namespace
