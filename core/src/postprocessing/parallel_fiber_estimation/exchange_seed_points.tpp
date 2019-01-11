#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
exchangeSeedPointsBeforeTracing(int nRanksZ, int rankZNo, bool streamlineDirectionUpwards, std::vector<Vec3> &seedPoints)
{
  // determine if previously set seedPoints are used or if they are received from neighbouring rank
  LOG(DEBUG) << "rankZNo: " << rankZNo << ", streamlineDirectionUpwards: " << streamlineDirectionUpwards;
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
    std::vector<double> receiveBuffer(seedPoints.size()*3);
    MPIUtility::handleReturnValue(MPI_Recv(receiveBuffer.data(), receiveBuffer.size(), MPI_DOUBLE, neighbourRankNo,
                                            0, currentRankSubset_->mpiCommunicator(), MPI_STATUS_IGNORE), "MPI_Recv");

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

    //LOG(DEBUG) << "send " << seedPoints.size() << " seed points to rank " << neighbourRankNo << ": " << seedPoints;

    // send seed points
    MPIUtility::handleReturnValue(MPI_Send(sendBuffer.data(), sendBuffer.size(), MPI_DOUBLE, neighbourRankNo,
                                            0, currentRankSubset_->mpiCommunicator()), "MPI_Send");
  }

}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
exchangeSeedPointsAfterTracing(int nRanksZ, int rankZNo, bool streamlineDirectionUpwards, std::vector<Vec3> &seedPoints, std::vector<std::vector<Vec3>> &streamlinePoints)
{
  if (nRanksZ == 1)
    return;

  // send end points of streamlines to next rank that continues the streamline
  if (nRanksZ > 1 && rankZNo != nRanksZ-1 && rankZNo != 0)
  {
    // fill send buffer
    std::vector<double> sendBuffer(seedPoints.size()*3);
    for (int streamlineIndex = 0; streamlineIndex < seedPoints.size(); streamlineIndex++)
    {
      for (int i = 0; i < 3; i++)
      {
        sendBuffer[streamlineIndex*3+i] = streamlinePoints[streamlineIndex].back()[i];
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

}

} // namespace
