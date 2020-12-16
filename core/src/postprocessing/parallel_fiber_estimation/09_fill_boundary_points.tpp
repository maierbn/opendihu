#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
fillBoundaryPoints(std::array<std::vector<std::vector<Vec3>>,4> &boundaryPoints, std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &boundaryPointsSubdomain,
                 std::array<std::vector<Vec3>,4> &cornerStreamlines, std::array<std::array<std::vector<bool>,4>,8> &boundaryPointsSubdomainAreValid,
                 std::array<bool,4> &subdomainIsAtBoundary)
{
  LOG(DEBUG) << "fillBoundaryPoints";
  // this method fills in the points that are on the outer boundary of the domain ("muscle surface")

  int nRanksZ = meshPartition_->nRanks(2);
  int rankZNo = meshPartition_->ownRankPartitioningIndex(2);

  double zRangeTotal = topZClip_ - bottomZClip_;
  double zRangeCurrentLevel = zRangeTotal / nRanksZ;
  double bottomZClip = bottomZClip_ + zRangeCurrentLevel*rankZNo;
  double topZClip = bottomZClip_ + zRangeCurrentLevel*(rankZNo+1);

  std::array<std::pair<int,int>,4> startEndCornerIndices = {
    std::pair<int,int>{0,2}, std::pair<int,int>{1,3}, std::pair<int,int>{0,1}, std::pair<int,int>{2,3}
  };
  //                                      0- 0+ 1- 1+
  std::array<int,4> leftNeighbourFace = {Mesh::face_t::face1Plus, Mesh::face_t::face1Minus, Mesh::face_t::face0Minus, Mesh::face_t::face0Plus};
  std::array<int,4> rightNeighbourFace = {Mesh::face_t::face1Minus, Mesh::face_t::face1Plus, Mesh::face_t::face0Plus, Mesh::face_t::face0Minus};

  // loop over z levels
  double currentZ;
  for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZNew_; zLevelIndex++)
  {
    currentZ = bottomZClip + double(zLevelIndex) / (nBoundaryPointsZNew_-1) * (topZClip - bottomZClip);

    // communication variables and buffers
    std::array<int,4> leftNeighbourRankNo;
    std::array<int,4> rightNeighbourRankNo;

    std::array<Vec3,4> leftBoundaryPoint;
    std::array<Vec3,4> leftForeignBoundaryPoint;
    std::array<Vec3,4> rightBoundaryPoint;
    std::array<Vec3,4> rightForeignBoundaryPoint;

    MPI_Request request;
    std::vector<MPI_Request> mpiRequests;

    //   ^ --(1+)-> ^
    // ^ 0-         0+
    // | | --(1-)-> |
    // +-->

    // determine start and end points for each face that is at the boundary

    for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
    {
      // if "face" is at the outer muscle boundary
      if (subdomainIsAtBoundary[face])
      {
        Vec3 startPoint, endPoint;

        // get start and end point of horizontal line on boundary that is to be found, in counterclockwise direction
        // by default, use the previous boundary point of the old boundary
        if (zLevelIndex % 2 == 0)
        {
          startPoint = boundaryPoints[face][zLevelIndex/2][0];
          endPoint = boundaryPoints[face][zLevelIndex/2][nBoundaryPointsXNew_-1];
        }
        else
        {
          startPoint = 0.5 * (boundaryPoints[face][zLevelIndex/2][0] + boundaryPoints[face][zLevelIndex/2+1][0]);
          endPoint = 0.5 * (boundaryPoints[face][zLevelIndex/2][nBoundaryPointsXNew_-1] + boundaryPoints[face][zLevelIndex/2+1][nBoundaryPointsXNew_-1]);
        }

        // if available, use traced streamlines along corners [0],[1],[2],[3] to determine startPoint and endPoint for boundary mesh
        //  [2]           [3]
        //    ^ --(1+)-> ^
        // ^  0-         0+
        // |  | --(1-)-> |
        // |[0]           [1]
        // +-->

        int beginIndex = startEndCornerIndices[face].first;
        int endIndex = startEndCornerIndices[face].second;

        // if the respective corner streamlines were corretly traced, use them instead of the old boundary points
        if (cornerStreamlines[beginIndex].size() == nBoundaryPointsZNew_)
        {
          LOG(DEBUG) << "cornerStreamline " << beginIndex << " is valid, use new start point " << cornerStreamlines[beginIndex][zLevelIndex] << " instead of " << startPoint;
          startPoint = cornerStreamlines[beginIndex][zLevelIndex];
        }

        if (cornerStreamlines[endIndex].size() == nBoundaryPointsZNew_)
        {
          LOG(DEBUG) << "cornerStreamline " << endIndex << " is valid, use new end point " << cornerStreamlines[endIndex][zLevelIndex] << " instead of " << endPoint;
          endPoint = cornerStreamlines[endIndex][zLevelIndex];
        }

        // assure that the neighbour rank uses the same boundary point
        // this is done by taking the average point with the neighbouring rank
        leftBoundaryPoint[face] = startPoint;
        rightBoundaryPoint[face] = endPoint;

        // left and right refers to the beginning and end of the arrow, which is left/right or bottom/top of this scheme:
        //    ^ --(1+)-> ^
        // ^  0-         0+
        // |  | --(1-)-> |
        // +-->
        if (face == (int)Mesh::face_t::face1Plus || face == (int)Mesh::face_t::face0Minus)
        {
          leftBoundaryPoint[face] = endPoint;
          rightBoundaryPoint[face] = startPoint;
        }

        leftForeignBoundaryPoint[face] = leftBoundaryPoint[face];
        rightForeignBoundaryPoint[face] = rightBoundaryPoint[face];

        leftNeighbourRankNo[face] = meshPartition_->neighbourRank((Mesh::face_or_edge_t)leftNeighbourFace[face]);
        if (leftNeighbourRankNo[face] != -1)
        {
          LOG(DEBUG) << "zLevelIndex " << zLevelIndex << ", face " << Mesh::getString((Mesh::face_t)face)
            << ", send boundaryPoint " << leftBoundaryPoint[face] << " to left rank " << leftNeighbourRankNo[face] << ", receive from left rank";

          // send and receive left boundary point
          MPIUtility::handleReturnValue(MPI_Isend(leftBoundaryPoint[face].data(), 3, MPI_DOUBLE, leftNeighbourRankNo[face], 1,
                                                  currentRankSubset_->mpiCommunicator(), &request), "MPI_Isend");
          mpiRequests.push_back(request);

          MPIUtility::handleReturnValue(MPI_Irecv(leftForeignBoundaryPoint[face].data(), 3, MPI_DOUBLE, leftNeighbourRankNo[face], 1,
                                                  currentRankSubset_->mpiCommunicator(), &request), "MPI_Irecv");
          mpiRequests.push_back(request);
        }
        else
        {
          LOG(DEBUG) << "zLevelIndex " << zLevelIndex << ", face " << Mesh::getString((Mesh::face_t)face) << ", has no left neighbour";
        }

        rightNeighbourRankNo[face] = meshPartition_->neighbourRank((Mesh::face_or_edge_t)rightNeighbourFace[face]);
        if (rightNeighbourRankNo[face] != -1)
        {
          LOG(DEBUG) << "zLevelIndex " << zLevelIndex << ", face " << Mesh::getString((Mesh::face_t)face)
            << ", send boundaryPoint " << rightBoundaryPoint[face] << " to right rank " << rightNeighbourRankNo[face] << ", receive from right rank";

          // send and receive right boundary point
          MPIUtility::handleReturnValue(MPI_Isend(rightBoundaryPoint[face].data(), 3, MPI_DOUBLE, rightNeighbourRankNo[face], 1,
                                                  currentRankSubset_->mpiCommunicator(), &request), "MPI_Isend");
          mpiRequests.push_back(request);

          MPIUtility::handleReturnValue(MPI_Irecv(rightForeignBoundaryPoint[face].data(), 3, MPI_DOUBLE, rightNeighbourRankNo[face], 1,
                                                  currentRankSubset_->mpiCommunicator(), &request), "MPI_Irecv");
          mpiRequests.push_back(request);
        }
        else
        {
          LOG(DEBUG) << "zLevelIndex " << zLevelIndex << ", face " << Mesh::getString((Mesh::face_t)face) << ", has no right neighbour";
        }
      }
    }

    // wait for communication to finish
    if (!mpiRequests.empty())
      MPIUtility::handleReturnValue(MPI_Waitall(mpiRequests.size(), mpiRequests.data(), MPI_STATUSES_IGNORE), "MPI_Waitall");

    for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
    {
      // if "face" is at the outer muscle boundary
      if (subdomainIsAtBoundary[face])
      {
        // start end end points are determined here for the face
        //  [2]           [3]
        //    ^ --(1+)-> ^
        // ^  0-         0+
        // |  | --(1-)-> |
        // |[0]           [1]
        // +-->

        LOG(DEBUG) << "left: neighbourRankNo: " << leftNeighbourRankNo[face]
          << ", boundaryPoint: " << leftBoundaryPoint[face] << ", foreignBoundaryPoint: " << leftForeignBoundaryPoint[face];
        LOG(DEBUG) << "right: neighbourRankNo: " << rightNeighbourRankNo[face]
          << ", boundaryPoint: " << rightBoundaryPoint[face] << ", foreignBoundaryPoint: " << rightForeignBoundaryPoint[face];

        // compute boundary points as average of the two corner streamlines of the adjacent processes
        leftBoundaryPoint[face] = 0.5*(leftBoundaryPoint[face] + leftForeignBoundaryPoint[face]);
        rightBoundaryPoint[face] = 0.5*(rightBoundaryPoint[face] + rightForeignBoundaryPoint[face]);

        LOG(DEBUG) << "new leftBoundaryPoint: " << leftBoundaryPoint << ", new rightBoundaryPoint: " << rightBoundaryPoint;

        Vec3 startPoint, endPoint;
        if (face == (int)Mesh::face_t::face1Plus || face == (int)Mesh::face_t::face0Minus)
        {
          endPoint = leftBoundaryPoint[face];
          startPoint = rightBoundaryPoint[face];
        }
        else
        {
          startPoint = leftBoundaryPoint[face];
          endPoint = rightBoundaryPoint[face];
        }
        // here, the points startPoint and endPoint define the line segment on the boundary that has to be determined now

#ifndef NDEBUG
#ifdef STL_OUTPUT
#ifdef STL_OUTPUT_VERBOSE
        if (zLevelIndex == 0 || (zLevelIndex >= nBoundaryPointsZ_-1 && zLevelIndex <= nBoundaryPointsZ_+1))
        {
          std::vector<Vec3> p;
          p.push_back(startPoint);
          p.push_back(endPoint);
          std::stringstream s;
          s << "07_start_end_point_face_" << Mesh::getString((Mesh::face_t)face) << "_z" << zLevelIndex;
          PyObject_CallFunction(functionOutputPoints_, "s i i O f", s.str().c_str(), currentRankSubset_->ownRankNo(), level_,
                                PythonUtility::convertToPython<std::vector<Vec3>>::get(p), 0.5);
          PythonUtility::checkForError();
        }
#endif
#endif
#endif
        LOG(DEBUG) << "z: " << zLevelIndex << ", face " << Mesh::getString(Mesh::face_t(face))
          << ", startPoint: " << startPoint << ", endPoint: " << endPoint << ", currentZ: " << currentZ
          << ", nBoundaryPointsXNew_: " << nBoundaryPointsXNew_;

        // determine boundary points between startPoint and endPoint, with same zLevel value
        // call stl_create_rings.create_ring_section_mesh

        PyObject *startPointPy = PythonUtility::convertToPython<Vec3>::get(startPoint);
        PyObject *endPointPy = PythonUtility::convertToPython<Vec3>::get(endPoint);
        PyObject *loopSectionPy = PyObject_CallFunction(functionCreateRingSection_, "s O O f i", inputMeshFilename_.c_str(), startPointPy, endPointPy, currentZ, nBoundaryPointsXNew_);
        PythonUtility::checkForError();
        //  create_ring_section(input_filename, start_point, end_point, z_value, n_points)
        assert(loopSectionPy);

        std::vector<Vec3> loopSection = PythonUtility::convertFromPython<std::vector<Vec3>>::get(loopSectionPy);
        if (loopSection.size() != nBoundaryPointsXNew_)
          LOG(ERROR) << "Ring section from python script contains only " << loopSection.size() << " points, " << nBoundaryPointsXNew_ << " requested.";
        assert(loopSection.size() == nBoundaryPointsXNew_);

        //LOG(DEBUG) << "got result " << loopSection;

        //       1+
        //    _________
        //   *    x    *
        //   |  2 | 3  |
        //0- |x___|___x|  0+
        //   |    |    |
        //   |  0 | 1  |
        //   *____x____*
        //
        //       1-
        // x = points that are not set by streamlines, but need to be set by the boundary

        // determine affected subdomains
        int subdomainIndex0 = 0;
        int subdomainIndex1 = 0;
        int facePerpendicular0 = 0;
        int facePerpendicular1 = 0;
        int pointIndexEdge = 0;

        if (face == (int)Mesh::face_t::face0Minus)
        {
          subdomainIndex0 = 0;
          subdomainIndex1 = 2;
          facePerpendicular0 = (int)Mesh::face_t::face1Plus;
          facePerpendicular1 = (int)Mesh::face_t::face1Minus;
          pointIndexEdge = 0;
        }
        else if (face == (int)Mesh::face_t::face0Plus)
        {
          subdomainIndex0 = 1;
          subdomainIndex1 = 3;
          facePerpendicular0 = (int)Mesh::face_t::face1Plus;
          facePerpendicular1 = (int)Mesh::face_t::face1Minus;
          pointIndexEdge = nBoundaryPointsX_-1;
        }
        else if (face == (int)Mesh::face_t::face1Minus)
        {
          subdomainIndex0 = 0;
          subdomainIndex1 = 1;
          facePerpendicular0 = (int)Mesh::face_t::face0Plus;
          facePerpendicular1 = (int)Mesh::face_t::face0Minus;
          pointIndexEdge = 0;
        }
        else if (face == (int)Mesh::face_t::face1Plus)
        {
          subdomainIndex0 = 2;
          subdomainIndex1 = 3;
          facePerpendicular0 = (int)Mesh::face_t::face0Plus;
          facePerpendicular1 = (int)Mesh::face_t::face0Minus;
          pointIndexEdge = nBoundaryPointsX_-1;
        }

        // bottom subdomain
        if (zLevelIndex <= nBoundaryPointsZ_-1)
        {
          int zLevelIndexSubdomain = zLevelIndex;
          assert(zLevelIndexSubdomain < nBoundaryPointsZ_);
          //LOG(DEBUG) << "save for subdomains " << subdomainIndex0 << " and " << subdomainIndex1 << ", zLevelIndex " << zLevelIndexSubdomain;
          int i = 0;
          for (std::vector<Vec3>::const_iterator loopPointIter = loopSection.begin(); loopPointIter != loopSection.begin()+nBoundaryPointsX_; loopPointIter++, i++)
          {

            assert(i < nBoundaryPointsX_);
            boundaryPointsSubdomain[subdomainIndex0][face][zLevelIndexSubdomain][i] = *loopPointIter;
          }

          i = 0;
          for (std::vector<Vec3>::const_iterator loopPointIter = loopSection.begin()+nBoundaryPointsX_-1; loopPointIter != loopSection.end(); loopPointIter++, i++)
          {
            assert(i < nBoundaryPointsX_);
            boundaryPointsSubdomain[subdomainIndex1][face][zLevelIndexSubdomain][i] = *loopPointIter;
          }

          // set end point of line of inner cross that does not have a traced streamline there, because it is at the boundary ("x" in figure)
          boundaryPointsSubdomain[subdomainIndex0][facePerpendicular0][zLevelIndexSubdomain][pointIndexEdge] = loopSection[nBoundaryPointsX_-1];
          boundaryPointsSubdomain[subdomainIndex1][facePerpendicular1][zLevelIndexSubdomain][pointIndexEdge] = loopSection[nBoundaryPointsX_-1];

          boundaryPointsSubdomainAreValid[subdomainIndex0][facePerpendicular0][pointIndexEdge] = true;
          boundaryPointsSubdomainAreValid[subdomainIndex1][facePerpendicular1][pointIndexEdge] = true;

          // set neighbouring boundary point at the end ("*" in figure)
          boundaryPointsSubdomain[subdomainIndex0][facePerpendicular1][zLevelIndexSubdomain][pointIndexEdge] = loopSection[0];
          boundaryPointsSubdomain[subdomainIndex1][facePerpendicular0][zLevelIndexSubdomain][pointIndexEdge] = loopSection[nBoundaryPointsXNew_-1];

          boundaryPointsSubdomainAreValid[subdomainIndex0][facePerpendicular1][pointIndexEdge] = true;
          boundaryPointsSubdomainAreValid[subdomainIndex1][facePerpendicular0][pointIndexEdge] = true;

#ifndef NDEBUG
#if 0
#ifdef STL_OUTPUT
          s.str("");
          s << "08_boundary_points_filled_face_" << Mesh::getString((Mesh::face_t)face) << "_z_" << zLevelIndexSubdomain << "_subdomain_" << subdomainIndex0;
          PyObject_CallFunction(functionOutputPoints_, "s i i O f", s.str().c_str(), currentRankSubset_->ownRankNo(), level_,
                                PythonUtility::convertToPython<std::vector<Vec3>>::get(boundaryPointsSubdomain[subdomainIndex0][face][zLevelIndexSubdomain]), 0.2);
          PythonUtility::checkForError();
          s.str("");
          s << "08_boundary_points_filled_face_" << Mesh::getString((Mesh::face_t)face) << "_z_" << zLevelIndexSubdomain << "_subdomain_" << subdomainIndex1;
          PyObject_CallFunction(functionOutputPoints_, "s i i O f", s.str().c_str(), currentRankSubset_->ownRankNo(), level_,
                                PythonUtility::convertToPython<std::vector<Vec3>>::get(boundaryPointsSubdomain[subdomainIndex1][face][zLevelIndexSubdomain]), 0.2);
          PythonUtility::checkForError();
#endif
#endif
#endif
          //std::copy(loopSection.begin(), loopSection.begin()+nBoundaryPointsX_, boundaryPointsSubdomain[subdomainIndex0][face][zLevelIndexSubdomain].begin());
          //std::copy(loopSection.begin()+nBoundaryPointsX_-1, loopSection.end(), boundaryPointsSubdomain[subdomainIndex1][face][zLevelIndexSubdomain].begin());
        }

        // top subdomain
        if (zLevelIndex >= nBoundaryPointsZ_-1)
        {
          subdomainIndex0 += 4;
          subdomainIndex1 += 4;
          int zLevelIndexSubdomain = zLevelIndex - (nBoundaryPointsZ_-1);
          assert(zLevelIndexSubdomain < nBoundaryPointsZ_);

          //LOG(DEBUG) << "save for subdomains " << subdomainIndex0 << " and " << subdomainIndex1 << ", zLevelIndex " << zLevelIndexSubdomain;
          //std::copy(loopSection.begin(), loopSection.begin()+nBoundaryPointsX_, boundaryPointsSubdomain[subdomainIndex0][face][zLevelIndexSubdomain].begin());
          //std::copy(loopSection.begin()+nBoundaryPointsX_-1, loopSection.end(), boundaryPointsSubdomain[subdomainIndex1][face][zLevelIndexSubdomain].begin());

          int i = 0;
          for (std::vector<Vec3>::const_iterator loopPointIter = loopSection.begin(); loopPointIter != loopSection.begin()+nBoundaryPointsX_; loopPointIter++, i++)
          {
            assert(i < nBoundaryPointsX_);
            boundaryPointsSubdomain[subdomainIndex0][face][zLevelIndexSubdomain][i] = *loopPointIter;
          }

          i = 0;
          for (std::vector<Vec3>::const_iterator loopPointIter = loopSection.begin()+nBoundaryPointsX_-1; loopPointIter != loopSection.end(); loopPointIter++, i++)
          {
            assert(i < nBoundaryPointsX_);
            boundaryPointsSubdomain[subdomainIndex1][face][zLevelIndexSubdomain][i] = *loopPointIter;
          }

          // set "center point" of loop section ("x" in figure)
          boundaryPointsSubdomain[subdomainIndex0][facePerpendicular0][zLevelIndexSubdomain][pointIndexEdge] = loopSection[nBoundaryPointsX_-1];
          boundaryPointsSubdomain[subdomainIndex1][facePerpendicular1][zLevelIndexSubdomain][pointIndexEdge] = loopSection[nBoundaryPointsX_-1];

          boundaryPointsSubdomainAreValid[subdomainIndex0][facePerpendicular0][pointIndexEdge] = true;
          boundaryPointsSubdomainAreValid[subdomainIndex1][facePerpendicular1][pointIndexEdge] = true;

          // set neighbouring boundary point at the end ("*" in figure)
          boundaryPointsSubdomain[subdomainIndex0][facePerpendicular1][zLevelIndexSubdomain][pointIndexEdge] = loopSection[0];
          boundaryPointsSubdomain[subdomainIndex1][facePerpendicular0][zLevelIndexSubdomain][pointIndexEdge] = loopSection[nBoundaryPointsXNew_-1];

          boundaryPointsSubdomainAreValid[subdomainIndex0][facePerpendicular1][pointIndexEdge] = true;
          boundaryPointsSubdomainAreValid[subdomainIndex1][facePerpendicular0][pointIndexEdge] = true;
#ifndef NDEBUG
#if 0
#ifdef STL_OUTPUT
          s.str("");
          s << "08_boundary_points_filled_face_" << Mesh::getString((Mesh::face_t)face) << "_z_" << zLevelIndexSubdomain << "_subdomain_" << subdomainIndex0;
          PyObject_CallFunction(functionOutputPoints_, "s i i O f", s.str().c_str(), currentRankSubset_->ownRankNo(), level_,
                                PythonUtility::convertToPython<std::vector<Vec3>>::get(boundaryPointsSubdomain[subdomainIndex0][face][zLevelIndexSubdomain]), 0.2);
          PythonUtility::checkForError();
          s.str("");
          s << "08_boundary_points_filled_face_" << Mesh::getString((Mesh::face_t)face) << "_z_" << zLevelIndexSubdomain << "_subdomain_" << subdomainIndex1;
          PyObject_CallFunction(functionOutputPoints_, "s i i O f", s.str().c_str(), currentRankSubset_->ownRankNo(), level_,
                                PythonUtility::convertToPython<std::vector<Vec3>>::get(boundaryPointsSubdomain[subdomainIndex1][face][zLevelIndexSubdomain]), 0.2);
          PythonUtility::checkForError();
#endif
#endif
#endif
        }

      }  // if subdomain is at boundary
    }  // for face
  }  // for z level

#ifndef NDEBUG
#ifdef STL_OUTPUT
#ifdef STL_OUTPUT_VERBOSE
  for (int subdomainIndex = 0; subdomainIndex < 8; subdomainIndex++)
  {
    std::stringstream s;
    s << "08_boundary_points_filled_subdomain_" << subdomainIndex;
    PyObject_CallFunction(functionOutputBoundaryPoints_, "s i i O f", s.str().c_str(), currentRankSubset_->ownRankNo(), level_,
                          PythonUtility::convertToPython<std::array<std::vector<std::vector<Vec3>>,4>>::get(boundaryPointsSubdomain[subdomainIndex]), 0.1);
    PythonUtility::checkForError();
  }
#endif
#endif
#endif

  // output points for debugging
#if 0
#ifndef NDEBUG
  std::ofstream file("points.csv", std::ios::out | std::ios::trunc);
  assert (file.is_open());
  file << "subdomainIndex;faceNo;zLevelIndex;;points" << std::endl;
  for (int subdomainIndex = 0; subdomainIndex < 8; subdomainIndex++)
  {
    for (int faceNo = (int)Mesh::face_t::face0Minus; faceNo <= (int)Mesh::face_t::face1Plus; faceNo++)
    {
      for (int zLevelIndex = 0; zLevelIndex < boundaryPointsSubdomain[subdomainIndex][faceNo].size(); zLevelIndex++)
      {
        file << subdomainIndex << ";" << faceNo << ";" << zLevelIndex << ";;";
        for (int pointIndex = 0; pointIndex < boundaryPointsSubdomain[subdomainIndex][faceNo][zLevelIndex].size(); pointIndex++)
        {
          file << boundaryPointsSubdomain[subdomainIndex][faceNo][zLevelIndex][pointIndex][0] << ";" << boundaryPointsSubdomain[subdomainIndex][faceNo][zLevelIndex][pointIndex][1] << ";" << boundaryPointsSubdomain[subdomainIndex][faceNo][zLevelIndex][pointIndex][2] << ";";
        }
        file << std::endl;
      }
    }
  }
  file.close();
#endif
#endif
}

} // namespace
