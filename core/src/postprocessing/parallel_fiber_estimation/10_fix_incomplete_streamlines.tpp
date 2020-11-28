#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
fixIncompleteStreamlines(std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &boundaryPointsSubdomain,
                         std::array<std::array<std::vector<bool>,4>,8> &boundaryPointsSubdomainAreValid, bool streamlineDirectionUpwards,
                         const std::array<bool,4> &subdomainIsAtBoundary, std::array<std::vector<std::vector<Vec3>>,4> boundaryPoints)
{
  // std::array<std::vector<std::vector<Vec3>>,4> boundaryPoints;    // [face_t][z-level][pointIndex]
  // boundaryPointsSubdomain[subdomainIndex][face][zLevelIndex][pointIndex]
  // boundaryPointsSubdomainAreValid[subdomainIndex][face][pointIndex]

#ifndef NDEBUG
  LOG(DEBUG) << "valid streamlines:";
  // output which streamlines are valid
  for (int subdomainIndex = 0; subdomainIndex < 8; subdomainIndex++)
  {
    for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
    {
      for (int pointIndex = 0; pointIndex < nBoundaryPointsX_; pointIndex++)
      {
        std::stringstream str;

        std::vector<Vec3> points;
        for (int i = 0; i < nBoundaryPointsZ_; i++)
        {
          points.push_back(boundaryPointsSubdomain[subdomainIndex][face][i][pointIndex]);
        }
        str << points;

        LOG(DEBUG) << "subdomain " << subdomainIndex << ", face " << Mesh::getString((Mesh::face_t)face) << ", streamline " << pointIndex << " valid: "
          << std::boolalpha << boundaryPointsSubdomainAreValid[subdomainIndex][face][pointIndex] << ", points: " << str.str();
      }
    }
  }
#endif

  // if there are streamlines at the edge between two processes' subdomains that are valid on one process and invalid on the other, send them from the valid process to the invalid
  communicateEdgeStreamlines(boundaryPointsSubdomain, boundaryPointsSubdomainAreValid);

  // fill invalid streamlines at corners from boundary points
  fixStreamlinesCorner(boundaryPointsSubdomain, boundaryPointsSubdomainAreValid, subdomainIsAtBoundary, boundaryPoints);

  // fill invalid streamlines, loop over the bottom 4 subdomains, the top are considered at the same iteration
  fixStreamlinesInterior(boundaryPointsSubdomain, boundaryPointsSubdomainAreValid, streamlineDirectionUpwards);

#ifndef NDEBUG
  LOG(DEBUG) << "after fixing, valid streamlines: ";

  int nInvalid = 0;
  // output which streamlines are valid
  for (int subdomainIndex = 0; subdomainIndex < 8; subdomainIndex++)
  {
    for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
    {
      for (int pointIndex = 0; pointIndex < nBoundaryPointsX_; pointIndex++)
      {

        if (!boundaryPointsSubdomainAreValid[subdomainIndex][face][pointIndex])
        {
          nInvalid++;
        }
        LOG(DEBUG) << "subdomain " << subdomainIndex << ", face " << Mesh::getString((Mesh::face_t)face) << ", streamline " << pointIndex << " valid: "
          << std::boolalpha << boundaryPointsSubdomainAreValid[subdomainIndex][face][pointIndex];
      }
    }
  }
  LOG(DEBUG) << "n invalid: " << nInvalid;

  if (nInvalid > 0)
    LOG(ERROR) << "There are " << nInvalid << " invalid streamlines! Level " << level_;
#endif
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
communicateEdgeStreamlines(std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &boundaryPointsSubdomain,
                           std::array<std::array<std::vector<bool>,4>,8> &boundaryPointsSubdomainAreValid)
{
  // if there are streamlines at the edge between two processes' subdomains that are valid on one process and invalid on the other, send them from the valid process to the invalid
  LOG(DEBUG) << "communicateEdgeStreamlines";
  //MPI_Barrier(currentRankSubset_->mpiCommunicator());

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
    int neighbourRankNo = meshPartition_->neighbourRank((Mesh::face_or_edge_t)face);

    // if there is no neighbour in that direction
    if (neighbourRankNo == -1)
      continue;

    sendBufferValidity[face].resize(nBoundaryPointsXNew_, 0);
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
      int pointIndexEnd = nBoundaryPointsXNew_;
      switch ((Mesh::face_t)face)
      {
        case Mesh::face_t::face0Minus:
        case Mesh::face_t::face0Plus:
          pointIndexStart = subdomainY*(nBoundaryPointsX_-1);
          pointIndexEnd = (subdomainY+1)*(nBoundaryPointsX_-1)+1;
          break;

        case Mesh::face_t::face1Minus:
        case Mesh::face_t::face1Plus:
          pointIndexStart = subdomainX*(nBoundaryPointsX_-1);
          pointIndexEnd = (subdomainX+1)*(nBoundaryPointsX_-1)+1;
          break;
        default:
          break;
      };

      LOG(DEBUG) << "face " << Mesh::getString((Mesh::face_t)face) << " subdomain " << subdomainIndex << ", streamlines [" << pointIndexStart << "," << pointIndexEnd << ") nBoundaryPointsXNew_: " << nBoundaryPointsXNew_;

      for (int pointIndex = pointIndexStart; pointIndex != pointIndexEnd; pointIndex++)
      {
        LOG(DEBUG) << "   " << pointIndex << ": " << std::boolalpha << boundaryPointsSubdomainAreValid[subdomainIndex][face][pointIndex-pointIndexStart];
        sendBufferValidity[face][pointIndex] = boundaryPointsSubdomainAreValid[subdomainIndex][face][pointIndex-pointIndexStart];
      }
    }

#ifndef NDEBUG
    std::stringstream str;
    for (int pointIndex = 0; pointIndex < nBoundaryPointsXNew_; pointIndex++)
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
    receiveBufferValidity[face].resize(nBoundaryPointsXNew_, 0);
    MPIUtility::handleReturnValue(MPI_Irecv(receiveBufferValidity[face].data(), nBoundaryPointsXNew_, MPI_CHAR,
                                            neighbourRankNo, tag, currentRankSubset_->mpiCommunicator(), &receiveRequest), "MPI_Irecv");
    receiveRequests.push_back(receiveRequest);

    // post non-blocking send call
    tag = this->currentRankSubset_->ownRankNo()*10000 + neighbourRankNo*100 + level_ + 8;
    LOG(DEBUG) << "send validity to rank " << neighbourRankNo << ", (" << tag << ")";
    MPI_Request sendRequest;
    MPIUtility::handleReturnValue(MPI_Isend(sendBufferValidity[face].data(), nBoundaryPointsXNew_, MPI_CHAR,
                                            neighbourRankNo, tag, currentRankSubset_->mpiCommunicator(), &sendRequest), "MPI_Isend");
    sendRequests.push_back(sendRequest);
  }

  // wait for non-blocking communication to finish
  if (!sendRequests.empty())
    MPIUtility::handleReturnValue(MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");
  if (!receiveRequests.empty())
    MPIUtility::handleReturnValue(MPI_Waitall(receiveRequests.size(), receiveRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

  LOG(DEBUG) << "waitall (" << sendRequests.size() << " send requests, " << receiveRequests.size() << " receiveRequests) complete";

  for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
  {
    int neighbourRankNo = meshPartition_->neighbourRank((Mesh::face_or_edge_t)face);

    // if there is no neighbour in that direction
    if (neighbourRankNo == -1)
      continue;

    std::stringstream str;
    for (int pointIndex = 0; pointIndex < nBoundaryPointsXNew_; pointIndex++)
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
  std::stringstream logMessage;

  int nStreamlinesFixed = 0;
  // loop over faces
  for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
  {
    int neighbourRankNo = meshPartition_->neighbourRank((Mesh::face_or_edge_t)face);

    // if there is no neighbour in that direction
    if (neighbourRankNo == -1)
      continue;

    sendBufferStreamline[face].resize(nBoundaryPointsXNew_);
    receiveBufferStreamline[face].resize(nBoundaryPointsXNew_);

    for (int pointIndex = 0; pointIndex != nBoundaryPointsXNew_; pointIndex++)
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
      if (pointIndex >= nBoundaryPointsX_-1)
      {
        subdomainIndex = subdomainIndex1;
        streamlineIndex = pointIndex - (nBoundaryPointsX_-1);
      }

      if (sendBufferValidity[face][pointIndex] && !receiveBufferValidity[face][pointIndex])
      {
        nStreamlinesFixed++;

        logMessage << "    face " << Mesh::getString((Mesh::face_t)face) << " subdomain " << subdomainIndex << ", streamline " << streamlineIndex
          << ", valid on own subdomain (rank " << currentRankSubset_->ownRankNo() << "), invalid on neighbor subdomain (rank " << neighbourRankNo << ")" << std::endl;

        // send streamline to neighbouring process
        sendBufferStreamline[face][pointIndex].resize(nBoundaryPointsZNew_*3);

        // bottom subdomain
        for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZ_; zLevelIndex++)
        {
          for (int i = 0; i < 3; i++)
          {
            sendBufferStreamline[face][pointIndex][3*zLevelIndex + i] = boundaryPointsSubdomain[subdomainIndex][face][zLevelIndex][streamlineIndex][i];
          }
        }

        // top subdomain
        for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZ_; zLevelIndex++)
        {
          for (int i = 0; i < 3; i++)
          {
            sendBufferStreamline[face][pointIndex][3*(zLevelIndex + (nBoundaryPointsZ_-1)) + i] = boundaryPointsSubdomain[subdomainIndex+4][face][zLevelIndex][streamlineIndex][i];
          }
        }

        // post non-blocking send call
        int tag = this->currentRankSubset_->ownRankNo()*10000 + neighbourRankNo*100 + level_ + 9;
        LOG(DEBUG) << "send streamline to rank " << neighbourRankNo << ", (" << tag << ")";
        MPI_Request sendRequest;
        assert(sendBufferStreamline[face][pointIndex].size() == nBoundaryPointsZNew_*3);
        MPIUtility::handleReturnValue(MPI_Isend(sendBufferStreamline[face][pointIndex].data(), nBoundaryPointsZNew_*3, MPI_DOUBLE,
                                                neighbourRankNo, tag, currentRankSubset_->mpiCommunicator(), &sendRequest), "MPI_Isend");
        sendRequests.push_back(sendRequest);
      }
      else if (!sendBufferValidity[face][pointIndex] && receiveBufferValidity[face][pointIndex])
      {
        // receive streamline from neighbouring process
        receiveBufferStreamline[face][pointIndex].resize(nBoundaryPointsZNew_*3);

        // post non-blocking receive call
        int tag = this->currentRankSubset_->ownRankNo()*100 + neighbourRankNo*10000 + level_ + 9;
        LOG(DEBUG) << "receive streamline from rank " << neighbourRankNo << ", (" << tag << ")";
        MPI_Request receiveRequest;
        assert (sizeof(bool) == 1);
        MPIUtility::handleReturnValue(MPI_Irecv(receiveBufferStreamline[face][pointIndex].data(), nBoundaryPointsZNew_*3, MPI_DOUBLE,
                                                neighbourRankNo, tag, currentRankSubset_->mpiCommunicator(), &receiveRequest), "MPI_Irecv");
        receiveRequests.push_back(receiveRequest);
      }
    }
  }

  // wait for non-blocking communication to finish
  if (!sendRequests.empty())
    MPIUtility::handleReturnValue(MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");
  if (!receiveRequests.empty())
    MPIUtility::handleReturnValue(MPI_Waitall(receiveRequests.size(), receiveRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

  LOG(DEBUG) << "waitall (" << sendRequests.size() << " send requests, " << receiveRequests.size() << " receiveRequests) complete";

  // assign received streamline data
  // loop over faces
  for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
  {
    int neighbourRankNo = meshPartition_->neighbourRank((Mesh::face_or_edge_t)face);

    // if there is no neighbour in that direction
    if (neighbourRankNo == -1)
      continue;

    for (int pointIndex = 0; pointIndex != nBoundaryPointsXNew_; pointIndex++)
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
        if (pointIndex >= nBoundaryPointsX_-1)
        {
          subdomainIndex = subdomainIndex1;
          streamlineIndex = pointIndex - (nBoundaryPointsX_-1);
        }

        // bottom subdomain
        for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZ_; zLevelIndex++)
        {
          Vec3 point;
          for (int i = 0; i < 3; i++)
          {
            point[i] = receiveBufferStreamline[face][pointIndex][3*zLevelIndex + i];
          }
          boundaryPointsSubdomain[subdomainIndex][face][zLevelIndex][streamlineIndex] = point;
        }

        // top subdomain
        for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZ_; zLevelIndex++)
        {
          Vec3 point;
          for (int i = 0; i < 3; i++)
          {
            point[i] = receiveBufferStreamline[face][pointIndex][3*(zLevelIndex + (nBoundaryPointsZ_-1)) + i];
          }
          boundaryPointsSubdomain[subdomainIndex+4][face][zLevelIndex][streamlineIndex] = point;
        }

        // set fixed streamlines to valid
        boundaryPointsSubdomainAreValid[subdomainIndex][face][streamlineIndex] = true;
        boundaryPointsSubdomainAreValid[subdomainIndex+4][face][streamlineIndex] = true;
      }
    }
  }

  // log invalid key fibers to log file
  if (nStreamlinesFixed > 0)
  {
    std::ofstream file;
    std::string logFilename = "out/log_fixed_streamlines.txt";
    OutputWriter::Generic::openFile(file, logFilename, true);
    file << currentRankSubset_->ownRankNo() << ": l=" << level_ << " communicateEdgeStreamlines, nFixed: " << nStreamlinesFixed << "\n"
      << logMessage.str();
    file.close();
  }
  LOG(DEBUG) << "end of communicateEdgeStreamlines";
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
fixStreamlinesCorner(std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &boundaryPointsSubdomain,
                     std::array<std::array<std::vector<bool>,4>,8> &boundaryPointsSubdomainAreValid,
                     const std::array<bool,4> &subdomainIsAtBoundary, std::array<std::vector<std::vector<Vec3>>,4> boundaryPoints)
{
  // fill invalid streamlines at corners from boundary points
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
    std::tuple<int,int,int>{1, Mesh::face_t::face1Minus, nBoundaryPointsX_-1}, std::tuple<int,int,int>{1, Mesh::face_t::face0Plus, 0},
    std::tuple<int,int,int>{2, Mesh::face_t::face0Minus, nBoundaryPointsX_-1}, std::tuple<int,int,int>{2, Mesh::face_t::face1Plus, 0},
    std::tuple<int,int,int>{3, Mesh::face_t::face1Plus, nBoundaryPointsX_-1}, std::tuple<int,int,int>{3, Mesh::face_t::face0Plus, nBoundaryPointsX_-1}
  };
  LOG(DEBUG) << "fixStreamlinesCorner";
  int nStreamlinesFixed = 0;
  std::stringstream logMessage;

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
    if (!boundaryPointsSubdomainAreValid[subdomainIndex0][face0][pointIndex0]
      || !boundaryPointsSubdomainAreValid[subdomainIndex1][face1][pointIndex1])
    {
      LOG(ERROR) << "found invalid streamline at corner, "
        << "subdomain " << subdomainIndex0 << " " << Mesh::getString((Mesh::face_t)face0) << "," << Mesh::getString((Mesh::face_t)face1);
      fixCorner = true;
    }

    // if it is an interior boundary definitely fix it, because it may be corrupted by sending it between processes
    /*if (!subdomainIsAtBoundary[face0] && !subdomainIsAtBoundary[face1])
    {
      LOG(ERROR) << "found streamline at corner in interior, "
        << "subdomain " << subdomainIndex0 << " " << Mesh::getString((Mesh::face_t)face0) << "," << Mesh::getString((Mesh::face_t)face1);
      fixCorner = true;
    }*/

    // if at least one of the two instances of the corner streamline is invalid
    if (fixCorner)
    {
      nStreamlinesFixed++;
      // derive the corner streamline from the initial boundary points

      // debugging output for log file
      logMessage << "    corner streamline of new subdomain " << subdomainIndex0 << ", face " << Mesh::getString((Mesh::face_t)face0)
        << "," << Mesh::getString((Mesh::face_t)face1) << ", derive from boundary points" << std::endl;

      // debugging output
      LOG(DEBUG) << " subdomain " << subdomainIndex0 << " face " << Mesh::getString((Mesh::face_t)face0)
        << " pointIndex " << pointIndex0 << " valid: " << std::boolalpha << boundaryPointsSubdomainAreValid[subdomainIndex0][face0][pointIndex0];
      LOG(DEBUG) << "           " << subdomainIndex1 << " face " << Mesh::getString((Mesh::face_t)face1)
        << " pointIndex " << pointIndex1 << " valid: " << std::boolalpha << boundaryPointsSubdomainAreValid[subdomainIndex1][face1][pointIndex1];
      LOG(DEBUG) << "derive the corner streamline from the initial boundary points";

      // lower subdomain
      // loop over points of streamline from bottom to top
      for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZ_; zLevelIndex++)
      {
        int pointIndexBoundaryPoints = 0;
        if (pointIndex0 > 0)
          pointIndexBoundaryPoints = nBoundaryPointsXNew_ - 1;

        LOG(DEBUG) << " z for boundaryPointsSubdomain: " << zLevelIndex << ", pointIndexBoundaryPoints: " << pointIndexBoundaryPoints;

        // determine point from boundaryPoints
        Vec3 streamlinePoint;
        Vec3 streamlinePointUpperSubdomain;
        if (zLevelIndex % 2 == 0)
        {
          streamlinePoint = boundaryPoints[face0][zLevelIndex/2][pointIndexBoundaryPoints];

          int zLevelIndexUpperSubdomain = int(zLevelIndex/2)+(nBoundaryPointsZ_-1)/2;
          streamlinePointUpperSubdomain = boundaryPoints[face0][zLevelIndexUpperSubdomain][pointIndexBoundaryPoints];
        }
        else
        {
          streamlinePoint = 0.5 * (boundaryPoints[face0][int(zLevelIndex/2)][pointIndexBoundaryPoints]
            + boundaryPoints[face0][int(zLevelIndex/2)+1][pointIndexBoundaryPoints]);

          int zLevelIndexUpperSubdomain = int(zLevelIndex/2)+(nBoundaryPointsZ_-1)/2;
          streamlinePointUpperSubdomain =
            0.5 * (boundaryPoints[face0][zLevelIndexUpperSubdomain][pointIndexBoundaryPoints]
                 + boundaryPoints[face0][zLevelIndexUpperSubdomain+1][pointIndexBoundaryPoints]);
        }

        LOG(DEBUG) << "streamlinePoints: lower: " << streamlinePoint << ", upper: " << streamlinePointUpperSubdomain;

        // lower subdomain
        boundaryPointsSubdomain[subdomainIndex0][face0][zLevelIndex][pointIndex0] = streamlinePoint;
        boundaryPointsSubdomain[subdomainIndex1][face1][zLevelIndex][pointIndex1] = streamlinePoint;

        // upper subdomain
        boundaryPointsSubdomain[subdomainIndex0+4][face0][zLevelIndex][pointIndex0] = streamlinePointUpperSubdomain;
        boundaryPointsSubdomain[subdomainIndex1+4][face1][zLevelIndex][pointIndex1] = streamlinePointUpperSubdomain;
      }

      // set fixed streamlines to valid
      boundaryPointsSubdomainAreValid[subdomainIndex0][face0][pointIndex0] = true;
      boundaryPointsSubdomainAreValid[subdomainIndex1][face1][pointIndex1] = true;
      boundaryPointsSubdomainAreValid[subdomainIndex0+4][face0][pointIndex0] = true;
      boundaryPointsSubdomainAreValid[subdomainIndex1+4][face1][pointIndex1] = true;
    }
  }
  LOG(DEBUG) << "fixStreamlinesCorner: " << nStreamlinesFixed << " fixed";

  // log invalid key fibers to log file
  if (nStreamlinesFixed > 0)
  {
    std::ofstream file;
    std::string logFilename = "out/log_fixed_streamlines.txt";
    OutputWriter::Generic::openFile(file, logFilename, true);
    file << currentRankSubset_->ownRankNo() << ": l=" << level_ << " fixStreamlinesCorner, nFixed: " << nStreamlinesFixed << "\n"
      << logMessage.str();
    file.close();
  }
  
  LOG(DEBUG) << "end of fixStreamlinesCorner";
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
fixStreamlinesInterior(std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &boundaryPointsSubdomain,
                       std::array<std::array<std::vector<bool>,4>,8> &boundaryPointsSubdomainAreValid, bool streamlineDirectionUpwards)
{
  LOG(DEBUG) << "fixStreamlinesInterior";
  int nStreamlinesFixed = 0;
  std::stringstream logMessage;

  // fill invalid streamlines, loop over the bottom 4 subdomains, the top are considered at the same iteration
  for (int subdomainIndex = 0; subdomainIndex < 4; subdomainIndex++)
  {
    for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
    {
      int lastValid = -1;
      bool lastStreamlinesWereInvalid = false;
      for (int pointIndex = 0; pointIndex < nBoundaryPointsX_; pointIndex++)
      {
        // if current streamline is valid
        if (boundaryPointsSubdomainAreValid[subdomainIndex][face][pointIndex])
        {
          if (lastStreamlinesWereInvalid)
          {
            LOG(DEBUG) << "subdomain " << subdomainIndex << ", face " << Mesh::getString((Mesh::face_t)face)
              << ", invalid streamlines from " << lastValid+1 << " to " << pointIndex-1 << ", lastValid: " << lastValid << ", nextValid: " << pointIndex;

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
                seedPointZLevelIndex = nBoundaryPointsZ_-1;
                seedPointSubdomainIndex += 4;
              }

              LOG(DEBUG) << "invalidStreamline " << invalidStreamlineIndex << ", seedPointZLevelIndex: " << seedPointZLevelIndex
                << ", seedPointSubdomainIndex: " << seedPointSubdomainIndex;

              // if there is a previous valid streamline
              if (lastValid != -1)
              {
                nStreamlinesFixed++;
                LOG(DEBUG) << "number z levels subdomain: " << boundaryPointsSubdomain[subdomainIndex][face].size();
                LOG(DEBUG) << "number z levels subdomain with seed point: " << boundaryPointsSubdomain[seedPointSubdomainIndex][face].size();
                LOG(DEBUG) << "number streamlines on subdomain with seed point at seed point (i.e. number of seed points): "
                  << boundaryPointsSubdomain[seedPointSubdomainIndex][face][seedPointZLevelIndex].size();

                Vec3 seedPointInvalid = boundaryPointsSubdomain[subdomainIndex][face][0][invalidStreamlineIndex];
                Vec3 seedPointLastValid = boundaryPointsSubdomain[seedPointSubdomainIndex][face][seedPointZLevelIndex][lastValid];
                Vec3 seedPointCurrent = boundaryPointsSubdomain[seedPointSubdomainIndex][face][seedPointZLevelIndex][pointIndex];

                LOG(DEBUG) << "seedPointInvalid: " << seedPointInvalid;
                LOG(DEBUG) << "seedPointLastValid: " << seedPointLastValid;
                LOG(DEBUG) << "seedPointCurrent: " << seedPointCurrent;

                LOG(DEBUG) << "lastValid: " << lastValid << ", streamline: ";
                for (int i = 0; i < nBoundaryPointsZ_; i++)
                {
                  LOG(DEBUG) << "zIndex: " << i << ": " << boundaryPointsSubdomain[seedPointSubdomainIndex][face][i][lastValid];
                }

                LOG(DEBUG) << "(next valid) pointIndex: " << pointIndex << ", current streamline: ";
                for (int i = 0; i < nBoundaryPointsZ_; i++)
                {
                  LOG(DEBUG) << "zIndex: " << i << ": " << boundaryPointsSubdomain[seedPointSubdomainIndex][face][i][pointIndex];
                }

                double alpha = MathUtility::norm<3>(seedPointInvalid - seedPointLastValid) / MathUtility::norm<3>(seedPointCurrent - seedPointLastValid);

                LOG(DEBUG) << "alpha: " << alpha;

                // debugging output for log file
                logMessage << "    invalidStreamline " << invalidStreamlineIndex << ", seedPointZLevelIndex: " << seedPointZLevelIndex
                  << ", seedPointSubdomainIndex: " << seedPointSubdomainIndex << ", alpha: " << alpha;

                if (alpha < 0 || alpha > 1)
                  LOG(WARNING) << "Interpolating invalid streamline in subdomain " << subdomainIndex
                    << ", face " << Mesh::getString((Mesh::face_t)face) << " with alpha = " << alpha << ", alpha should be in [0,1].";

                //loop over points of streamline from bottom to top
                for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZ_; zLevelIndex++)
                {
                  boundaryPointsSubdomain[subdomainIndex][face][zLevelIndex][invalidStreamlineIndex]
                    = (1.-alpha) * boundaryPointsSubdomain[subdomainIndex][face][zLevelIndex][lastValid]
                      + alpha * boundaryPointsSubdomain[subdomainIndex][face][zLevelIndex][pointIndex];

                  boundaryPointsSubdomain[subdomainIndex+4][face][zLevelIndex][invalidStreamlineIndex]
                    = (1.-alpha) * boundaryPointsSubdomain[subdomainIndex+4][face][zLevelIndex][lastValid]
                      + alpha * boundaryPointsSubdomain[subdomainIndex+4][face][zLevelIndex][pointIndex];
                }
                LOG(DEBUG) << "interpolate in subdomain " << subdomainIndex << ", face " << Mesh::getString((Mesh::face_t)face)
                  << " streamline " << invalidStreamlineIndex << " from " << lastValid << " and " << pointIndex << ", alpha: " << alpha;

                // set fixed streamline to valid
                boundaryPointsSubdomainAreValid[subdomainIndex][face][invalidStreamlineIndex] = true;
                boundaryPointsSubdomainAreValid[subdomainIndex+4][face][invalidStreamlineIndex] = true;
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

  // log invalid key fibers to log file
  if (nStreamlinesFixed > 0)
  {
    std::ofstream file;
    std::string logFilename = "out/log_fixed_streamlines.txt";
    OutputWriter::Generic::openFile(file, logFilename, true);
    file << currentRankSubset_->ownRankNo() << ": l=" << level_ << " fixStreamlinesInterior, nFixed: " << nStreamlinesFixed << "\n"
      << logMessage.str();
    file.close();
  }
  
  LOG(DEBUG) << "end of fixStreamlinesInterior";
}

} // namespace
