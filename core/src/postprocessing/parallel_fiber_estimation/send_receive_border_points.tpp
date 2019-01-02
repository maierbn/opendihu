#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
sendBorderPoints(std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &borderPointsSubdomain, std::vector<MPI_Request> &sendRequests)
{
  sendRequests.resize(8);

  // send border points to subdomains
  for (int subdomainIndex = 0; subdomainIndex < 8; subdomainIndex++)
  {
    // fill send buffer
    std::vector<double> sendBuffer;
    for (int faceNo = (int)Mesh::face_t::face0Minus; faceNo <= (int)Mesh::face_t::face1Plus; faceNo++)
    {
      for (int zLevelIndex = 0; zLevelIndex < borderPointsSubdomain[subdomainIndex][faceNo].size(); zLevelIndex++)
      {
        for (int pointIndex = 0; pointIndex < borderPointsSubdomain[subdomainIndex][faceNo][zLevelIndex].size(); pointIndex++)
        {
          for (int i = 0; i < 3; i++)
          {
            sendBuffer.push_back(borderPointsSubdomain[subdomainIndex][faceNo][zLevelIndex][pointIndex][i]);
          }
        }
      }
    }

    // determine rank so to which the border points should be sent to
    std::array<int,3> oldRankIndex;

    for (int i = 0; i < 3; i++)
    {
      oldRankIndex[i] = meshPartition_->ownRankPartitioningIndex(i);
    }

    std::array<int,3> subdomainRankIndex;
    subdomainRankIndex[0] = oldRankIndex[0]*2 + subdomainIndex % 2;
    subdomainRankIndex[1] = oldRankIndex[1]*2 + int((subdomainIndex % 4) / 2);
    subdomainRankIndex[2] = oldRankIndex[2]*2 + int(subdomainIndex / 4);

    int subdomainRankNo = subdomainRankIndex[2]*(nRanksPerCoordinateDirection_[0]*nRanksPerCoordinateDirection_[1]) + subdomainRankIndex[1]*nRanksPerCoordinateDirection_[0] + subdomainRankIndex[0];

    // determine if the subdomain is at any border of the whole domain
    std::array<int,3> rankIndex;
    rankIndex[0] = subdomainIndex % nRanksPerCoordinateDirection_[0];
    rankIndex[1] = int((subdomainIndex % (nRanksPerCoordinateDirection_[0]*nRanksPerCoordinateDirection_[1])) / nRanksPerCoordinateDirection_[0]);
    rankIndex[2] = int(subdomainIndex / (nRanksPerCoordinateDirection_[0]*nRanksPerCoordinateDirection_[1]));

    std::array<bool,4> subdomainIsAtBorderNew;

    subdomainIsAtBorderNew[(int)Mesh::face_t::face0Minus] = rankIndex[0] == 0;
    subdomainIsAtBorderNew[(int)Mesh::face_t::face0Plus] = rankIndex[0] == nRanksPerCoordinateDirection_[0]-1;
    subdomainIsAtBorderNew[(int)Mesh::face_t::face1Minus] = rankIndex[1] == 0;
    subdomainIsAtBorderNew[(int)Mesh::face_t::face1Plus] = rankIndex[1] == nRanksPerCoordinateDirection_[1]-1;


    LOG(DEBUG) << "subdomain " << subdomainIndex << ", rankIndex: " << rankIndex << ", subdomainIsAtBorderNew: " << subdomainIsAtBorderNew << ", nRanksPerCoordinateDirection_: " << nRanksPerCoordinateDirection_;

    LOG(DEBUG) << "sendBuffer size: " << sendBuffer.size() << ", nBorderPointsX_:" << nBorderPointsX_ << ", nBorderPointsZ_:" << nBorderPointsZ_ << ", " << nBorderPointsX_*nBorderPointsZ_*3*4;

    // save subdomains to be send to a file
#ifdef WRITE_CHECKPOINT_BORDER_POINTS
    std::stringstream filename;
    filename << "checkpoint_borderPoints_subdomain_" << subdomainRankNo << ".csv";
    std::ofstream file(filename.str().c_str(), std::ios::out | std::ios::trunc);
    assert(file.is_open());

    // header
    file << borderPointsSubdomain[subdomainIndex][0].size() << ";" << borderPointsSubdomain[subdomainIndex][0][0].size() << ";";
    for (int i = 0; i < 4; i++)
    {
      file << (subdomainIsAtBorderNew[i]? "1" : "0") << ";";
    }
    file << std::endl;

    // data
    for (int faceNo = (int)Mesh::face_t::face0Minus; faceNo <= (int)Mesh::face_t::face1Plus; faceNo++)
    {
      for (int zLevelIndex = 0; zLevelIndex < borderPointsSubdomain[subdomainIndex][faceNo].size(); zLevelIndex++)
      {
        for (int pointIndex = 0; pointIndex < borderPointsSubdomain[subdomainIndex][faceNo][zLevelIndex].size(); pointIndex++)
        {
          for (int i = 0; i < 3; i++)
          {
            file << borderPointsSubdomain[subdomainRankNo][faceNo][zLevelIndex][pointIndex][i] << ";";
          }
        }
        file << std::endl;
      }
    }

    file.close();
    LOG(DEBUG) << " saved data for subdomain " << subdomainRankNo << " to file \"" << filename.str() << "\".";
#endif

    LOG(DEBUG) << "oldRankIndex: " << oldRankIndex << ", new subdomainIndex: " << subdomainIndex << ", subdomainRankIndex: " << subdomainRankIndex << ", send from rank " << currentRankSubset_->ownRankNo() << " to " << subdomainRankNo;
    MPIUtility::handleReturnValue(MPI_Isend(sendBuffer.data(), sendBuffer.size(), MPI_DOUBLE, subdomainRankNo, 0, currentRankSubset_->mpiCommunicator(), &sendRequests[subdomainIndex]), "MPI_Isend");

    LOG(DEBUG) << "subdomain " << subdomainIndex << " sendBuffer: " << sendBuffer;
  }
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
receiveBorderPoints(int nRanksPerCoordinateDirectionPreviously, std::array<std::vector<std::vector<Vec3>>,4> &borderPointsNew, std::array<bool,4> &subdomainIsAtBorderNew)
{
  int ownRankNo = currentRankSubset_->ownRankNo();

  // determine rank from which to receive the border points
  std::array<int,3> rankIndex;
  rankIndex[0] = ownRankNo % nRanksPerCoordinateDirection_[0];
  rankIndex[1] = int((ownRankNo % (nRanksPerCoordinateDirection_[0]*nRanksPerCoordinateDirection_[1])) / nRanksPerCoordinateDirection_[0]);
  rankIndex[2] = int(ownRankNo / (nRanksPerCoordinateDirection_[0]*nRanksPerCoordinateDirection_[1]));

  int rankToReceiveFrom = int(rankIndex[2]/2)*(nRanksPerCoordinateDirectionPreviously*nRanksPerCoordinateDirectionPreviously) + int(rankIndex[1]/2)*nRanksPerCoordinateDirectionPreviously + int(rankIndex[0]/2);

  int recvBufferSize = 4*nBorderPointsX_*nBorderPointsZ_*3;
  LOG(DEBUG) << "receive from rank " << rankToReceiveFrom << ", recvBufferSize: " << recvBufferSize << "(nBorderPointsX_: " << nBorderPointsX_ << ", nBorderPointsZ_: " << nBorderPointsZ_ << ")";

  std::vector<double> recvBuffer(recvBufferSize);
  MPIUtility::handleReturnValue(MPI_Recv(recvBuffer.data(), recvBuffer.size(), MPI_DOUBLE, rankToReceiveFrom, 0, currentRankSubset_->mpiCommunicator(), MPI_STATUS_IGNORE), "MPI_Recv");

  LOG(DEBUG) << "recv buffer: " << recvBuffer;

  // extract received values, assign them to borderPointsNew
  int recvBufferIndex = 0;
  for (int faceNo = (int)Mesh::face_t::face0Minus; faceNo <= (int)Mesh::face_t::face1Plus; faceNo++)
  {
    borderPointsNew[faceNo].resize(nBorderPointsZ_);
    for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
    {
      borderPointsNew[faceNo][zLevelIndex].resize(nBorderPointsX_);
      for (int pointIndex = 0; pointIndex < nBorderPointsX_; pointIndex++, recvBufferIndex+=3)
      {
        borderPointsNew[faceNo][zLevelIndex][pointIndex] = Vec3({recvBuffer[recvBufferIndex], recvBuffer[recvBufferIndex+1], recvBuffer[recvBufferIndex+2]});
      }
    }
  }
  assert(recvBufferIndex == recvBufferSize);
}

};  // namespace
