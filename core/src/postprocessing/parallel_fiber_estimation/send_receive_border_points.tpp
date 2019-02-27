#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
sendBorderPoints(std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &borderPointsSubdomain,
                 std::vector<std::vector<double>> &sendBuffers, std::vector<MPI_Request> &sendRequests)
{
  sendRequests.resize(8);
  sendBuffers.resize(8);

  int sendBufferSize = 4*nBorderPointsX_*nBorderPointsZ_*3;

  // send border points to subdomains
  for (int subdomainIndex = 0; subdomainIndex < 8; subdomainIndex++)
  {
    // fill send buffer
    int sendBufferIndex = 0;
    sendBuffers[subdomainIndex].resize(sendBufferSize);

    for (int faceNo = (int)Mesh::face_t::face0Minus; faceNo <= (int)Mesh::face_t::face1Plus; faceNo++)
    {
      for (int zLevelIndex = 0; zLevelIndex < borderPointsSubdomain[subdomainIndex][faceNo].size(); zLevelIndex++)
      {
        for (int pointIndex = 0; pointIndex < borderPointsSubdomain[subdomainIndex][faceNo][zLevelIndex].size(); pointIndex++)
        {
          for (int i = 0; i < 3; i++, sendBufferIndex++)
          {
            sendBuffers[subdomainIndex][sendBufferIndex] = borderPointsSubdomain[subdomainIndex][faceNo][zLevelIndex][pointIndex][i];
          }
        }
      }
    }

    assert(sendBufferIndex == sendBufferSize);

    // determine rank to which the border points should be sent to
    std::array<int,3> oldRankIndex;

    for (int i = 0; i < 3; i++)
    {
      oldRankIndex[i] = meshPartition_->ownRankPartitioningIndex(i);
    }
    
    LOG(DEBUG) << " --- ";
    LOG(DEBUG) << "subdomain " << subdomainIndex << ", oldRankIndex: " << oldRankIndex;

    std::array<int,3> subdomainRankIndex;
    subdomainRankIndex[0] = oldRankIndex[0]*2 + subdomainIndex % 2;
    subdomainRankIndex[1] = oldRankIndex[1]*2 + int((subdomainIndex % 4) / 2);
    subdomainRankIndex[2] = oldRankIndex[2]*2 + int(subdomainIndex / 4);

    int subdomainRankNo = subdomainRankIndex[2]*(nRanksPerCoordinateDirection_[0]*nRanksPerCoordinateDirection_[1])
      + subdomainRankIndex[1]*nRanksPerCoordinateDirection_[0] + subdomainRankIndex[0];

    LOG(DEBUG) << "subdomainRankIndex: " << subdomainRankIndex << ", subdomainRankNo: " << subdomainRankNo;

    // determine if the subdomain is at any border of the whole domain
    std::array<bool,4> subdomainIsAtBorderNew = {false,false,false,false};
    setSubdomainIsAtBorder(subdomainRankNo, subdomainIsAtBorderNew);

    LOG(DEBUG) << "subdomainIsAtBorderNew: " << subdomainIsAtBorderNew << ", sendBuffer size: " << sendBufferSize << ", nBorderPointsX_:" << nBorderPointsX_ << ", nBorderPointsZ_:"
      << nBorderPointsZ_ << ", " << nBorderPointsX_*nBorderPointsZ_*3*4;

    // save subdomains to be send to a file
#ifdef WRITE_CHECKPOINT_BORDER_POINTS
    std::stringstream filename;
    filename << "checkpoint_borderPoints_subdomain_" << subdomainRankNo << ".csv";
    LOG(DEBUG) << " filename \"" << filename.str() << "\".";
    std::ofstream file(filename.str().c_str(), std::ios::out | std::ios::trunc);
    if(!file.is_open())
    {
      LOG(WARNING) << "Could not open file \"" << filename.str() << "\".";
    }
    else 
    {
      LOG(DEBUG) << " opened file.";

      // header
      file << borderPointsSubdomain[subdomainIndex][0].size() << ";";
      LOG(DEBUG) << "a";
      if (borderPointsSubdomain[subdomainIndex][0].size() == 0)
      {
        file << "0;";
      }
      else
      {
        file << borderPointsSubdomain[subdomainIndex][0][0].size() << ";";
      }
      LOG(DEBUG) << "b";
      for (int i = 0; i < 4; i++)
      {
        file << (subdomainIsAtBorderNew[i]? "1" : "0") << ";";
        LOG(DEBUG) << "c" << i;
      }
      file << std::endl;
      LOG(DEBUG) << "d";
      
      // data
      for (int faceNo = (int)Mesh::face_t::face0Minus; faceNo <= (int)Mesh::face_t::face1Plus; faceNo++)
      {
        for (int zLevelIndex = 0; zLevelIndex < borderPointsSubdomain[subdomainIndex][faceNo].size(); zLevelIndex++)
        {
          for (int pointIndex = 0; pointIndex < borderPointsSubdomain[subdomainIndex][faceNo][zLevelIndex].size(); pointIndex++)
          {
            for (int i = 0; i < 3; i++)
            {
              LOG(DEBUG) << "e " << faceNo << " " << zLevelIndex << " " << pointIndex << " " << i;
              file << borderPointsSubdomain[subdomainRankNo][faceNo][zLevelIndex][pointIndex][i] << ";";
            }
          }
          file << std::endl;
        }
      }
      LOG(DEBUG) << "f";

      file.close();
      LOG(DEBUG) << " saved data for subdomain " << subdomainRankNo << " to file \"" << filename.str() << "\".";
    }
#endif

#ifdef FILE_COMMUNICATION
    std::stringstream filename;
    int ownRankNo = currentRankSubset_->ownRankNo();
    filename << "file_communication_borderPoints_from_" << ownRankNo << "_to_" << subdomainRankNo << ".txt";
    std::ofstream file(filename.str().c_str(), std::ios::out | std::ios::trunc);
    assert(file.is_open());

    // data
    for (int i = 0; i < sendBufferSize; i++)
    {
      file << sendBuffers[subdomainIndex][i] << " ";
    }

    file.close();
    LOG(DEBUG) << " file communication, saved to file \"" << filename.str() << "\".";
#endif

#ifndef FILE_COMMUNICATION
    int tag = currentRankSubset_->ownRankNo()*10000 + subdomainRankNo*100 + level_*10 + 4;
    LOG(DEBUG) << "oldRankIndex: " << oldRankIndex << ", new subdomainIndex: " << subdomainIndex << ", subdomainRankIndex: " << subdomainRankIndex
      << ", send from rank " << currentRankSubset_->ownRankNo() << " to " << subdomainRankNo << "(tag: " << tag << ")";
    MPIUtility::handleReturnValue(MPI_Isend(sendBuffers[subdomainIndex].data(), sendBufferSize, MPI_DOUBLE,
                                            subdomainRankNo, tag, currentRankSubset_->mpiCommunicator(), &sendRequests[subdomainIndex]), "MPI_Isend");
#endif
    LOG(DEBUG) << "subdomain " << subdomainIndex << " sendBuffer: " << sendBuffers[subdomainIndex];

#if 1
    std::stringstream debugFilename;
    debugFilename << "01_border_points_to_send_subdomain_" << subdomainIndex;
    PyObject_CallFunction(functionOutputBorderPoints_, "s i O f", debugFilename.str().c_str(), currentRankSubset_->ownRankNo(),
                          PythonUtility::convertToPython<std::array<std::vector<std::vector<Vec3>>,4>>::get(borderPointsSubdomain[subdomainIndex]), 0.1);
    PythonUtility::checkForError();
#endif

  }
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
receiveBorderPoints(int nRanksPerCoordinateDirectionPreviously, std::array<std::vector<std::vector<Vec3>>,4> &borderPointsNew,
                    std::array<bool,4> &subdomainIsAtBorderNew)
{
  int ownRankNo = currentRankSubset_->ownRankNo();

  // determine rank from which to receive the border points
  std::array<int,3> rankCoordinates;
  rankCoordinates[0] = ownRankNo % nRanksPerCoordinateDirection_[0];
  rankCoordinates[1] = int((ownRankNo % (nRanksPerCoordinateDirection_[0]*nRanksPerCoordinateDirection_[1])) / nRanksPerCoordinateDirection_[0]);
  rankCoordinates[2] = int(ownRankNo / (nRanksPerCoordinateDirection_[0]*nRanksPerCoordinateDirection_[1]));

  int rankToReceiveFrom = int(rankCoordinates[2]/2)*(nRanksPerCoordinateDirectionPreviously*nRanksPerCoordinateDirectionPreviously) + int(rankCoordinates[1]/2)*nRanksPerCoordinateDirectionPreviously + int(rankCoordinates[0]/2);

  int recvBufferSize = 4*nBorderPointsX_*nBorderPointsZ_*3;
  LOG(DEBUG) << "receive from rank " << rankToReceiveFrom << ", recvBufferSize: " << recvBufferSize
    << "(nBorderPointsX_: " << nBorderPointsX_ << ", nBorderPointsZ_: " << nBorderPointsZ_ << ")";

  std::vector<double> recvBuffer(recvBufferSize);
#ifdef FILE_COMMUNICATION
    std::stringstream filename;
    filename << "file_communication_borderPoints_from_" << rankToReceiveFrom << "_to_" << ownRankNo << ".txt";
    std::ifstream file;
    while(!file.is_open())
    {
      file.open(filename.str().c_str(), std::ios::in);
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }

    // data
    for (int i = 0; i < recvBuffer.size(); i++)
    {
      file >> recvBuffer[i];
    }

    file.close();
    LOG(DEBUG) << " file communication, loaded from file \"" << filename.str() << "\".";

#else
  MPI_Request receiveRequest;

  // receive border points
  int tag = currentRankSubset_->ownRankNo()*100 + rankToReceiveFrom*10000 + level_*10 + 4;
  MPIUtility::handleReturnValue(MPI_Irecv(recvBuffer.data(), recvBufferSize, MPI_DOUBLE,
                                          rankToReceiveFrom, tag, currentRankSubset_->mpiCommunicator(), &receiveRequest), "MPI_Irecv");
  MPIUtility::handleReturnValue(MPI_Wait(&receiveRequest, MPI_STATUS_IGNORE), "MPI_Wait");
  LOG(DEBUG) << " receive size " << recvBufferSize << ", tag: " << tag;

#endif
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

  // set subdomainIsAtBorderNew
  setSubdomainIsAtBorder(currentRankSubset_->ownRankNo(), subdomainIsAtBorderNew);

#if 1
    std::stringstream debugFilename;
    debugFilename << "02_border_points_received";
    PyObject_CallFunction(functionOutputBorderPoints_, "s i O f", debugFilename.str().c_str(), currentRankSubset_->ownRankNo(),
                          PythonUtility::convertToPython<std::array<std::vector<std::vector<Vec3>>,4>>::get(borderPointsNew), 0.1);
    PythonUtility::checkForError();
#endif

  assert(recvBufferIndex == recvBufferSize);
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
setSubdomainIsAtBorder(int rankNo, std::array<bool,4> &subdomainIsAtBorderNew)
{
  // determine if the subdomain is at which borders, from rank no
  std::array<int,3> rankCoordinates;
  rankCoordinates[0] = rankNo % nRanksPerCoordinateDirection_[0];
  rankCoordinates[1] = int((rankNo % (nRanksPerCoordinateDirection_[0]*nRanksPerCoordinateDirection_[1])) / nRanksPerCoordinateDirection_[0]);
  rankCoordinates[2] = int(rankNo / (nRanksPerCoordinateDirection_[0]*nRanksPerCoordinateDirection_[1]));

  subdomainIsAtBorderNew[(int)Mesh::face_t::face0Minus] = rankCoordinates[0] == 0;
  subdomainIsAtBorderNew[(int)Mesh::face_t::face0Plus] = rankCoordinates[0] == nRanksPerCoordinateDirection_[0]-1;
  subdomainIsAtBorderNew[(int)Mesh::face_t::face1Minus] = rankCoordinates[1] == 0;
  subdomainIsAtBorderNew[(int)Mesh::face_t::face1Plus] = rankCoordinates[1] == nRanksPerCoordinateDirection_[1]-1;

  LOG(DEBUG) << "subdomain rank " << rankNo << ", rankCoordinates: " << rankCoordinates << ", subdomainIsAtBorderNew: "
    << subdomainIsAtBorderNew << ", nRanksPerCoordinateDirection_: " << nRanksPerCoordinateDirection_;

}

} // namespace
