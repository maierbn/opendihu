#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
fixIncompleteStreamlines(std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &borderPointsSubdomain,
                         std::array<std::array<std::vector<bool>,4>,8> &borderPointsSubdomainAreValid, bool streamlineDirectionUpwards,
                         const std::array<bool,4> &subdomainIsAtBorder, std::array<std::vector<std::vector<Vec3>>,4> borderPoints)
{
  // std::array<std::vector<std::vector<Vec3>>,4> borderPoints;    // [face_t][z-level][pointIndex]
  // borderPointsSubdomain[subdomainIndex][face][zLevelIndex][pointIndex]
  // borderPointsSubdomainAreValid[subdomainIndex][face][pointIndex]

#ifndef NDEBUG
  LOG(DEBUG) << "valid streamlines:";
  // output which streamlines are valid
  for (int subdomainIndex = 0; subdomainIndex < 8; subdomainIndex++)
  {
    for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
    {
      for (int pointIndex = 0; pointIndex < nBorderPointsX_; pointIndex++)
      {
        std::stringstream str;

        std::vector<Vec3> points;
        for (int i = 0; i < nBorderPointsZ_; i++)
        {
          points.push_back(borderPointsSubdomain[subdomainIndex][face][i][pointIndex]);
        }
        str << points;

        LOG(DEBUG) << "subdomain " << subdomainIndex << ", face " << Mesh::getString((Mesh::face_t)face) << ", streamline " << pointIndex << " valid: "
          << std::boolalpha << borderPointsSubdomainAreValid[subdomainIndex][face][pointIndex] << ", points: " << str.str();
      }
    }
  }
#endif



  // if there are streamlines at the edge between two processes' subdomains that are valid on one process and invalid on the other, send them from the valid process to the invalid
  communicateEdgeStreamlines(borderPointsSubdomain, borderPointsSubdomainAreValid);

  // fill invalid streamlines at corners from border points
  fixStreamlinesCorner(borderPointsSubdomain, borderPointsSubdomainAreValid, subdomainIsAtBorder, borderPoints);

  // fill invalid streamlines, loop over the bottom 4 subdomains, the top are considered at the same iteration
  fixStreamlinesInterior(borderPointsSubdomain, borderPointsSubdomainAreValid, streamlineDirectionUpwards);

#ifndef NDEBUG
  LOG(DEBUG) << "after fixing, valid streamlines: ";

  int nInvalid = 0;
  // output which streamlines are valid
  for (int subdomainIndex = 0; subdomainIndex < 8; subdomainIndex++)
  {
    for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
    {
      for (int pointIndex = 0; pointIndex < nBorderPointsX_; pointIndex++)
      {

        if (!borderPointsSubdomainAreValid[subdomainIndex][face][pointIndex])
        {
          nInvalid++;
        }
        LOG(DEBUG) << "subdomain " << subdomainIndex << ", face " << Mesh::getString((Mesh::face_t)face) << ", streamline " << pointIndex << " valid: "
          << std::boolalpha << borderPointsSubdomainAreValid[subdomainIndex][face][pointIndex];
      }
    }
  }
  LOG(DEBUG) << "n invalid: " << nInvalid;
#endif
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
communicateEdgeStreamlines(std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &borderPointsSubdomain,
                           std::array<std::array<std::vector<bool>,4>,8> &borderPointsSubdomainAreValid)
{
  // if there are streamlines at the edge between two processes' subdomains that are valid on one process and invalid on the other, send them from the valid process to the invalid
  LOG(DEBUG) << "communicateEdgeStreamlines";
  MPI_Barrier(currentRankSubset_->mpiCommunicator());

  //
  //   ^ --(1+)-> ^   ^ --(1+)-> ^
  //   0-   [2]   0+  0-   [3]   0+
  //   | --(1-)-> |   | --(1-)-> |
  //
  //   ^ --(1+)-> ^   ^ --(1+)-> ^
  // ^ 0-   [0]   0+  0-   [1]   0+
  // | | --(1-)-> |   | --(1-)-> |
  // +-->

  std::vector<MPI_Request> receiveRequests;
  std::vector<MPI_Request> sendRequests;

  // communicate which streamlines are valid on the edges between processes
  std::array<std::vector<char>,4> sendBufferValidity;
  std::array<std::vector<char>,4> receiveBufferValidity;
  for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
  {
    int neighbourRankNo = meshPartition_->neighbourRank((Mesh::face_t)face);

    // if there is no neighbour in that direction
    if (neighbourRankNo == -1)
      continue;

    sendBufferValidity[face].resize(nBorderPointsXNew_, 0);
    for (int subdomainIndex = 0; subdomainIndex < 4; subdomainIndex++)
    {
      int subdomainX = subdomainIndex % 2;
      int subdomainY = int(subdomainIndex / 2);

      // only consider faces on the edge between processes' subdomains
      if (face == (int)Mesh::face_t::face0Minus && subdomainX != 0)
        continue;
      if (face == (int)Mesh::face_t::face0Plus && subdomainX != 1)
        continue;
      if (face == (int)Mesh::face_t::face1Minus && subdomainY != 0)
        continue;
      if (face == (int)Mesh::face_t::face1Plus && subdomainY != 1)
        continue;

      int pointIndexStart = 0;
      int pointIndexEnd = nBorderPointsXNew_;
      switch ((Mesh::face_t)face)
      {
        case Mesh::face_t::face0Minus:
        case Mesh::face_t::face0Plus:
          pointIndexStart = subdomainY*(nBorderPointsX_-1);
          pointIndexEnd = (subdomainY+1)*(nBorderPointsX_-1)+1;
          break;

        case Mesh::face_t::face1Minus:
        case Mesh::face_t::face1Plus:
          pointIndexStart = subdomainX*(nBorderPointsX_-1);
          pointIndexEnd = (subdomainX+1)*(nBorderPointsX_-1)+1;
          break;
        default:
          break;
      };

      LOG(DEBUG) << "face " << Mesh::getString((Mesh::face_t)face) << " subdomain " << subdomainIndex << ", streamlines [" << pointIndexStart << "," << pointIndexEnd << ") nBorderPointsXNew_: " << nBorderPointsXNew_;

      for (int pointIndex = pointIndexStart; pointIndex != pointIndexEnd; pointIndex++)
      {
        LOG(DEBUG) << "   " << pointIndex << ": " << std::boolalpha << borderPointsSubdomainAreValid[subdomainIndex][face][pointIndex-pointIndexStart];
        sendBufferValidity[face][pointIndex] = borderPointsSubdomainAreValid[subdomainIndex][face][pointIndex-pointIndexStart];
      }
    }

#ifndef NDEBUG
    std::stringstream str;
    for (int pointIndex = 0; pointIndex < nBorderPointsXNew_; pointIndex++)
    {
      if (pointIndex != 0)
        str << ",";
      str << (sendBufferValidity[face][pointIndex] == 0? "false" : "true");
    }
    LOG(DEBUG) << "sendBufferValidity " << Mesh::getString((Mesh::face_t)face) << ": [" << str.str() << "], neighbourRankNo: " << neighbourRankNo;
#endif

    // post non-blocking receive call
    int tag = this->currentRankSubset_->ownRankNo()*100 + neighbourRankNo*10000 + level_ + 8;
    LOG(DEBUG) << "receive validity from rank " << neighbourRankNo << ", (" << tag << ")";
    MPI_Request receiveRequest;
    receiveBufferValidity[face].resize(nBorderPointsXNew_, 0);
    MPIUtility::handleReturnValue(MPI_Irecv(receiveBufferValidity[face].data(), nBorderPointsXNew_, MPI_CHAR,
                                            neighbourRankNo, tag, currentRankSubset_->mpiCommunicator(), &receiveRequest), "MPI_Irecv");
    receiveRequests.push_back(receiveRequest);

    // post non-blocking send call
    tag = this->currentRankSubset_->ownRankNo()*10000 + neighbourRankNo*100 + level_ + 8;
    LOG(DEBUG) << "send validity to rank " << neighbourRankNo << ", (" << tag << ")";
    MPI_Request sendRequest;
    MPIUtility::handleReturnValue(MPI_Isend(sendBufferValidity[face].data(), nBorderPointsXNew_, MPI_CHAR,
                                            neighbourRankNo, tag, currentRankSubset_->mpiCommunicator(), &sendRequest), "MPI_Isend");
    sendRequests.push_back(sendRequest);
  }

  // wait for non-blocking communication to finish
  MPIUtility::handleReturnValue(MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");
  MPIUtility::handleReturnValue(MPI_Waitall(receiveRequests.size(), receiveRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

  LOG(DEBUG) << "waitall (" << sendRequests.size() << " send requests, " << receiveRequests.size() << " receiveRequests) complete";

  for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
  {
    int neighbourRankNo = meshPartition_->neighbourRank((Mesh::face_t)face);

    // if there is no neighbour in that direction
    if (neighbourRankNo == -1)
      continue;

    std::stringstream str;
    for (int pointIndex = 0; pointIndex < nBorderPointsXNew_; pointIndex++)
    {
      if (pointIndex != 0)
        str << ",";
      str << (receiveBufferValidity[face][pointIndex] == 0? "false" : "true");
    }
    LOG(DEBUG) << "receiveBufferValidity " << Mesh::getString((Mesh::face_t)face) << ": [" << str.str() << "], neighbourRankNo: " << neighbourRankNo;
  }

  // send and receive streamline data
  LOG(DEBUG) << "send and receive streamline data";
  sendRequests.clear();
  receiveRequests.clear();

  std::array<std::vector<std::vector<double>>,4> sendBufferStreamline;
  std::array<std::vector<std::vector<double>>,4> receiveBufferStreamline;

  int nStreamlinesFixed = 0;
  // loop over faces
  for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
  {
    int neighbourRankNo = meshPartition_->neighbourRank((Mesh::face_t)face);

    // if there is no neighbour in that direction
    if (neighbourRankNo == -1)
      continue;

    sendBufferStreamline[face].resize(nBorderPointsXNew_);
    receiveBufferStreamline[face].resize(nBorderPointsXNew_);

    for (int pointIndex = 0; pointIndex != nBorderPointsXNew_; pointIndex++)
    {
      int subdomainIndex0 = 0;
      int subdomainIndex1 = 0;
      switch ((Mesh::face_t)face)
      {
        case Mesh::face_t::face0Minus:
          subdomainIndex0 = 0;
          subdomainIndex1 = 2;
          break;
        case Mesh::face_t::face0Plus:
          subdomainIndex0 = 1;
          subdomainIndex1 = 3;
          break;
        case Mesh::face_t::face1Minus:
          subdomainIndex0 = 0;
          subdomainIndex1 = 1;
          break;
        case Mesh::face_t::face1Plus:
          subdomainIndex0 = 2;
          subdomainIndex1 = 3;
          break;
        default:
          break;
      };

      int subdomainIndex = subdomainIndex0;
      int streamlineIndex = pointIndex;
      if (pointIndex >= nBorderPointsX_-1)
      {
        subdomainIndex = subdomainIndex1;
        streamlineIndex = pointIndex - (nBorderPointsX_-1);
      }

      if (sendBufferValidity[face][pointIndex] && !receiveBufferValidity[face][pointIndex])
      {
        nStreamlinesFixed++;

        // send streamline to neighbouring process
        sendBufferStreamline[face][pointIndex].resize(nBorderPointsZNew_*3);

        // bottom subdomain
        for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
        {
          for (int i = 0; i < 3; i++)
          {
            sendBufferStreamline[face][pointIndex][3*zLevelIndex + i] = borderPointsSubdomain[subdomainIndex][face][zLevelIndex][streamlineIndex][i];
          }
        }

        // top subdomain
        for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
        {
          for (int i = 0; i < 3; i++)
          {
            sendBufferStreamline[face][pointIndex][3*(zLevelIndex + (nBorderPointsZ_-1)) + i] = borderPointsSubdomain[subdomainIndex+4][face][zLevelIndex][streamlineIndex][i];
          }
        }

        // post non-blocking send call
        int tag = this->currentRankSubset_->ownRankNo()*10000 + neighbourRankNo*100 + level_ + 9;
        LOG(DEBUG) << "send streamline to rank " << neighbourRankNo << ", (" << tag << ")";
        MPI_Request sendRequest;
        assert(sendBufferStreamline[face][pointIndex].size() == nBorderPointsZNew_*3);
        MPIUtility::handleReturnValue(MPI_Isend(sendBufferStreamline[face][pointIndex].data(), nBorderPointsZNew_*3, MPI_DOUBLE,
                                                neighbourRankNo, tag, currentRankSubset_->mpiCommunicator(), &sendRequest), "MPI_Isend");
        sendRequests.push_back(sendRequest);
      }
      else if (!sendBufferValidity[face][pointIndex] && receiveBufferValidity[face][pointIndex])
      {
        // receive streamline from neighbouring process
        receiveBufferStreamline[face][pointIndex].resize(nBorderPointsZNew_*3);

        // post non-blocking receive call
        int tag = this->currentRankSubset_->ownRankNo()*100 + neighbourRankNo*10000 + level_ + 9;
        LOG(DEBUG) << "receive streamline from rank " << neighbourRankNo << ", (" << tag << ")";
        MPI_Request receiveRequest;
        assert (sizeof(bool) == 1);
        MPIUtility::handleReturnValue(MPI_Irecv(receiveBufferStreamline[face][pointIndex].data(), nBorderPointsZNew_*3, MPI_DOUBLE,
                                                neighbourRankNo, tag, currentRankSubset_->mpiCommunicator(), &receiveRequest), "MPI_Irecv");
        receiveRequests.push_back(receiveRequest);
      }
    }
  }

  // wait for non-blocking communication to finish
  MPIUtility::handleReturnValue(MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");
  MPIUtility::handleReturnValue(MPI_Waitall(receiveRequests.size(), receiveRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

  LOG(DEBUG) << "waitall (" << sendRequests.size() << " send requests, " << receiveRequests.size() << " receiveRequests) complete";

  // assign received streamline data
  // loop over faces
  for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
  {
    int neighbourRankNo = meshPartition_->neighbourRank((Mesh::face_t)face);

    // if there is no neighbour in that direction
    if (neighbourRankNo == -1)
      continue;

    for (int pointIndex = 0; pointIndex != nBorderPointsXNew_; pointIndex++)
    {
      if (!sendBufferValidity[face][pointIndex] && receiveBufferValidity[face][pointIndex])
      {

        int subdomainIndex0 = 0;
        int subdomainIndex1 = 0;
        switch ((Mesh::face_t)face)
        {
          case Mesh::face_t::face0Minus:
            subdomainIndex0 = 0;
            subdomainIndex1 = 2;
            break;
          case Mesh::face_t::face0Plus:
            subdomainIndex0 = 1;
            subdomainIndex1 = 3;
            break;
          case Mesh::face_t::face1Minus:
            subdomainIndex0 = 0;
            subdomainIndex1 = 1;
            break;
          case Mesh::face_t::face1Plus:
            subdomainIndex0 = 2;
            subdomainIndex1 = 3;
            break;
          default:
            break;
        };

        int subdomainIndex = subdomainIndex0;
        int streamlineIndex = pointIndex;
        if (pointIndex >= nBorderPointsX_-1)
        {
          subdomainIndex = subdomainIndex1;
          streamlineIndex = pointIndex - (nBorderPointsX_-1);
        }

        // bottom subdomain
        for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
        {
          Vec3 point;
          for (int i = 0; i < 3; i++)
          {
            point[i] = receiveBufferStreamline[face][pointIndex][3*zLevelIndex + i];
          }
          borderPointsSubdomain[subdomainIndex][face][zLevelIndex][streamlineIndex] = point;
        }

        // top subdomain
        for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
        {
          Vec3 point;
          for (int i = 0; i < 3; i++)
          {
            point[i] = receiveBufferStreamline[face][pointIndex][3*(zLevelIndex + (nBorderPointsZ_-1)) + i];
          }
          borderPointsSubdomain[subdomainIndex+4][face][zLevelIndex][streamlineIndex] = point;
        }

        // set fixed streamlines to valid
        borderPointsSubdomainAreValid[subdomainIndex][face][streamlineIndex] = true;
        borderPointsSubdomainAreValid[subdomainIndex+4][face][streamlineIndex] = true;
      }
    }
  }
  LOG(DEBUG) << "communicateEdgeStreamlines: " << nStreamlinesFixed << " fixed";
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
fixStreamlinesCorner(std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &borderPointsSubdomain,
                     std::array<std::array<std::vector<bool>,4>,8> &borderPointsSubdomainAreValid,
                     const std::array<bool,4> &subdomainIsAtBorder, std::array<std::vector<std::vector<Vec3>>,4> borderPoints)
{
  // fill invalid streamlines at corners from border points
  //
  //   ^ --(1+)-> ^   ^ --(1+)-> ^
  //   0-   [2]   0+  0-   [3]   0+
  //   | --(1-)-> |   | --(1-)-> |
  //
  //   ^ --(1+)-> ^   ^ --(1+)-> ^
  // ^ 0-   [0]   0+  0-   [1]   0+
  // | | --(1-)-> |   | --(1-)-> |
  // +-->

  std::vector<std::tuple<int,int,int>> subdomainFacePointIndex = {
    std::tuple<int,int,int>{0, Mesh::face_t::face0Minus, 0}, std::tuple<int,int,int>{0, Mesh::face_t::face1Minus, 0},
    std::tuple<int,int,int>{1, Mesh::face_t::face1Minus, nBorderPointsX_-1}, std::tuple<int,int,int>{1, Mesh::face_t::face0Plus, 0},
    std::tuple<int,int,int>{2, Mesh::face_t::face0Minus, nBorderPointsX_-1}, std::tuple<int,int,int>{2, Mesh::face_t::face1Plus, 0},
    std::tuple<int,int,int>{3, Mesh::face_t::face1Plus, nBorderPointsX_-1}, std::tuple<int,int,int>{3, Mesh::face_t::face0Plus, nBorderPointsX_-1}
  };
  LOG(DEBUG) << "fixStreamlinesCorner";
  int nStreamlinesFixed = 0;

  // TODO: center point (not needed apparently)

  for (int i = 0; i < 8; i+=2)
  {
    int subdomainIndex0 = std::get<0>(subdomainFacePointIndex[i]);
    int face0 = std::get<1>(subdomainFacePointIndex[i]);
    int pointIndex0 = std::get<2>(subdomainFacePointIndex[i]);

    int subdomainIndex1 = std::get<0>(subdomainFacePointIndex[i+1]);
    int face1 = std::get<1>(subdomainFacePointIndex[i+1]);
    int pointIndex1 = std::get<2>(subdomainFacePointIndex[i+1]);

    bool fixCorner = false;

    // if at least one of the two instances of the corner streamline is invalid
    if (!borderPointsSubdomainAreValid[subdomainIndex0][face0][pointIndex0]
      || !borderPointsSubdomainAreValid[subdomainIndex1][face1][pointIndex1])
    {
      LOG(DEBUG) << "found invalid streamline at corner";
      fixCorner = true;
    }

    // if it is an interior border definitely fix it, because it may be corrupted by sending it between processes
    if (!subdomainIsAtBorder[face0] && !subdomainIsAtBorder[face1])
    {
      LOG(DEBUG) << "found streamline at corner in interior";
      fixCorner = true;
    }

    // if at least one of the two instances of the corner streamline is invalid
    if (fixCorner)
    {
      nStreamlinesFixed++;
      // derive the corner streamline from the initial border points

      LOG(DEBUG) << " subdomain " << subdomainIndex0 << " face " << Mesh::getString((Mesh::face_t)face0)
        << " pointIndex " << pointIndex0 << " valid: " << std::boolalpha << borderPointsSubdomainAreValid[subdomainIndex0][face0][pointIndex0];
      LOG(DEBUG) << "           " << subdomainIndex1 << " face " << Mesh::getString((Mesh::face_t)face1)
        << " pointIndex " << pointIndex1 << " valid: " << std::boolalpha << borderPointsSubdomainAreValid[subdomainIndex1][face1][pointIndex1];
      LOG(DEBUG) << "derive the corner streamline from the initial border points";

      // lower subdomain
      // loop over points of streamline from bottom to top
      for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
      {
        int pointIndexBorderPoints = 0;
        if (pointIndex0 > 0)
          pointIndexBorderPoints = nBorderPointsXNew_ - 1;

        LOG(DEBUG) << " z for borderPointsSubdomain: " << zLevelIndex << ", pointIndexBorderPoints: " << pointIndexBorderPoints;

        // determine point from borderPoints
        Vec3 streamlinePoint;
        Vec3 streamlinePointUpperSubdomain;
        if (zLevelIndex % 2 == 0)
        {
          streamlinePoint = borderPoints[face0][zLevelIndex/2][pointIndexBorderPoints];

          int zLevelIndexUpperSubdomain = int(zLevelIndex/2)+(nBorderPointsZ_-1)/2;
          streamlinePointUpperSubdomain = borderPoints[face0][zLevelIndexUpperSubdomain][pointIndexBorderPoints];
        }
        else
        {
          streamlinePoint = 0.5 * (borderPoints[face0][int(zLevelIndex/2)][pointIndexBorderPoints]
            + borderPoints[face0][int(zLevelIndex/2)+1][pointIndexBorderPoints]);

          int zLevelIndexUpperSubdomain = int(zLevelIndex/2)+(nBorderPointsZ_-1)/2;
          streamlinePointUpperSubdomain =
            0.5 * (borderPoints[face0][zLevelIndexUpperSubdomain][pointIndexBorderPoints]
                 + borderPoints[face0][zLevelIndexUpperSubdomain+1][pointIndexBorderPoints]);
        }

        LOG(DEBUG) << "streamlinePoints: lower: " << streamlinePoint << ", upper: " << streamlinePointUpperSubdomain;

        // lower subdomain
        borderPointsSubdomain[subdomainIndex0][face0][zLevelIndex][pointIndex0] = streamlinePoint;
        borderPointsSubdomain[subdomainIndex1][face1][zLevelIndex][pointIndex1] = streamlinePoint;

        // upper subdomain
        borderPointsSubdomain[subdomainIndex0+4][face0][zLevelIndex][pointIndex0] = streamlinePointUpperSubdomain;
        borderPointsSubdomain[subdomainIndex1+4][face1][zLevelIndex][pointIndex1] = streamlinePointUpperSubdomain;
      }

      // set fixed streamlines to valid
      borderPointsSubdomainAreValid[subdomainIndex0][face0][pointIndex0] = true;
      borderPointsSubdomainAreValid[subdomainIndex1][face1][pointIndex1] = true;
      borderPointsSubdomainAreValid[subdomainIndex0+4][face0][pointIndex0] = true;
      borderPointsSubdomainAreValid[subdomainIndex1+4][face1][pointIndex1] = true;
    }
  }
  LOG(DEBUG) << "fixStreamlinesCorner: " << nStreamlinesFixed << " fixed";
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
fixStreamlinesInterior(std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &borderPointsSubdomain,
                       std::array<std::array<std::vector<bool>,4>,8> &borderPointsSubdomainAreValid, bool streamlineDirectionUpwards)
{
  LOG(DEBUG) << "fixStreamlinesInterior";
  int nStreamlinesFixed = 0;

  // fill invalid streamlines, loop over the bottom 4 subdomains, the top are considered at the same iteration
  for (int subdomainIndex = 0; subdomainIndex < 4; subdomainIndex++)
  {
    for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
    {
      int lastValid = -1;
      bool lastStreamlinesWereInvalid = false;
      for (int pointIndex = 0; pointIndex < nBorderPointsX_; pointIndex++)
      {
        // if current streamline is valid
        if (borderPointsSubdomainAreValid[subdomainIndex][face][pointIndex])
        {
          if (lastStreamlinesWereInvalid)
          {
            LOG(DEBUG) << "subdomain " << subdomainIndex << ", face " << Mesh::getString((Mesh::face_t)face)
              << ", invalid streamlines " << lastValid+1 << " - " << pointIndex-1 << ", lastValid: " << lastValid << ", nextValid: " << pointIndex;

            // loop over streamlines between lastValid and pointIndex, these are all invalid, interpolate them
            for (int invalidStreamlineIndex = lastValid+1; invalidStreamlineIndex < pointIndex; invalidStreamlineIndex++)
            {
              int seedPointZLevelIndex = 0;
              int seedPointSubdomainIndex = subdomainIndex;
              if (streamlineDirectionUpwards)
              {
                seedPointZLevelIndex = 0;
              }
              else
              {
                seedPointZLevelIndex = nBorderPointsZ_-1;
                seedPointSubdomainIndex += 4;
              }

              LOG(DEBUG) << "invalidStreamline " << invalidStreamlineIndex << ", seedPointZLevelIndex: " << seedPointZLevelIndex
                << ", seedPointSubdomainIndex: " << seedPointSubdomainIndex;

              // if there is a previous valid streamline
              if (lastValid != -1)
              {
                nStreamlinesFixed++;
                LOG(DEBUG) << borderPointsSubdomain[subdomainIndex][face].size();
                LOG(DEBUG) << borderPointsSubdomain[seedPointSubdomainIndex][face].size();
                LOG(DEBUG) << borderPointsSubdomain[seedPointSubdomainIndex][face][seedPointZLevelIndex].size();

                Vec3 seedPointInvalid = borderPointsSubdomain[subdomainIndex][face][0][invalidStreamlineIndex];
                Vec3 seedPointLastValid = borderPointsSubdomain[seedPointSubdomainIndex][face][seedPointZLevelIndex][lastValid];
                Vec3 seedPointCurrent = borderPointsSubdomain[seedPointSubdomainIndex][face][seedPointZLevelIndex][pointIndex];

                LOG(DEBUG) << "seedPointInvalid: " << seedPointInvalid;
                LOG(DEBUG) << "seedPointLastValid: " << seedPointLastValid;
                LOG(DEBUG) << "seedPointCurrent: " << seedPointCurrent;

                LOG(DEBUG) << "lastValid: " << lastValid << ", streamline: ";
                for (int i = 0; i < nBorderPointsZ_; i++)
                {
                  LOG(DEBUG) << "zIndex: " << i << ": " << borderPointsSubdomain[seedPointSubdomainIndex][face][i][lastValid];
                }

                LOG(DEBUG) << "pointIndex: " << pointIndex << ", current streamline: ";
                for (int i = 0; i < nBorderPointsZ_; i++)
                {
                  LOG(DEBUG) << "zIndex: " << i << ": " << borderPointsSubdomain[seedPointSubdomainIndex][face][i][pointIndex];
                }

                double alpha = MathUtility::norm<3>(seedPointInvalid - seedPointLastValid) / MathUtility::norm<3>(seedPointCurrent - seedPointLastValid);

                LOG(DEBUG) << "alpha: " << alpha;

                //loop over points of streamline from bottom to top
                for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
                {
                  borderPointsSubdomain[subdomainIndex][face][zLevelIndex][invalidStreamlineIndex]
                    = (1.-alpha) * borderPointsSubdomain[subdomainIndex][face][zLevelIndex][lastValid]
                      + alpha * borderPointsSubdomain[subdomainIndex][face][zLevelIndex][pointIndex];

                  borderPointsSubdomain[subdomainIndex+4][face][zLevelIndex][invalidStreamlineIndex]
                    = (1.-alpha) * borderPointsSubdomain[subdomainIndex+4][face][zLevelIndex][lastValid]
                      + alpha * borderPointsSubdomain[subdomainIndex+4][face][zLevelIndex][pointIndex];
                }
                LOG(DEBUG) << "interpolate in subdomain " << subdomainIndex << ", face " << Mesh::getString((Mesh::face_t)face)
                  << " streamline " << invalidStreamlineIndex << " from " << lastValid << " and " << pointIndex << ", alpha: " << alpha;

                // set fixed streamline to valid
                borderPointsSubdomainAreValid[subdomainIndex][face][invalidStreamlineIndex] = true;
                borderPointsSubdomainAreValid[subdomainIndex+4][face][invalidStreamlineIndex] = true;
              }
              else
              {
                LOG(WARNING) << "Could not fix incomplete streamline on subdomain " << subdomainIndex
                  << ", face " << Mesh::getString((Mesh::face_t)face) << ", no " << invalidStreamlineIndex;
              }
            }
          }
          lastStreamlinesWereInvalid = false;
          lastValid = pointIndex;
        }
        else
        {
          lastStreamlinesWereInvalid = true;
        }
      }
    }
  }
  LOG(DEBUG) << "fixStreamlinesInterior: " << nStreamlinesFixed << " fixed";
}

} // namespace
