#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
fillBorderPoints(std::array<std::vector<std::vector<Vec3>>,4> &borderPoints, std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &borderPointsSubdomain,
                 std::array<bool,4> &subdomainIsAtBorder)
{
  PyObject *stlMeshPy = PyObject_CallFunction(functionGetStlMesh_, "s", stlFilename_.c_str());
  PythonUtility::checkForError();
  assert(stlMeshPy);

  int nRanksZ = meshPartition_->nRanks(2);
  int rankZNo = meshPartition_->ownRankPartitioningIndex(2);

  double zRangeTotal = topZClip_ - bottomZClip_;
  double zRangeCurrentLevel = zRangeTotal / nRanksZ;
  double bottomZClip = bottomZClip_ + zRangeCurrentLevel*rankZNo;
  double topZClip = bottomZClip_ + zRangeCurrentLevel*(rankZNo+1);

  // fill in points that are on the border of the domain
  // loop over z levels
  double currentZ;
  for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZNew_; zLevelIndex++)
  {
    currentZ = bottomZClip + double(zLevelIndex) / (nBorderPointsZNew_-1) * (topZClip - bottomZClip);

    Vec3 startPoint;
    Vec3 endPoint;

    //   ^ --(1+)-> ^
    // ^ 0-         0+
    // | | --(1-)-> |
    // +-->

    for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
    {
      if (subdomainIsAtBorder[face])
      {
        // get start and end point of horizontal line on boundary that is to be found
        if (zLevelIndex % 2 == 0)
        {
          startPoint = borderPoints[face][zLevelIndex/2][0];
          endPoint = borderPoints[face][zLevelIndex/2][nBorderPointsXNew_-1];
        }
        else
        {
          startPoint = 0.5 * (borderPoints[face][zLevelIndex/2][0] + borderPoints[face][zLevelIndex/2+1][0]);
          endPoint = 0.5 * (borderPoints[face][zLevelIndex/2][nBorderPointsXNew_-1] + borderPoints[face][zLevelIndex/2+1][nBorderPointsXNew_-1]);
        }


#ifndef NDEBUG
#ifdef STL_OUTPUT
#ifdef STL_OUTPUT_VERBOSE
        if (zLevelIndex == 0 || (zLevelIndex >= nBorderPointsZ_-1 && zLevelIndex <= nBorderPointsZ_+1))
        {
          std::vector<Vec3> p;
          p.push_back(startPoint);
          p.push_back(endPoint);
          std::stringstream s;
          s << "07_start_end_point_face_" << Mesh::getString((Mesh::face_t)face) << "_z" << zLevelIndex;
          PyObject_CallFunction(functionOutputPoints_, "s i O f", s.str().c_str(), currentRankSubset_->ownRankNo(),
                                PythonUtility::convertToPython<std::vector<Vec3>>::get(p), 0.5);
          PythonUtility::checkForError();
        }
#endif
#endif
#endif

        LOG(DEBUG) << "z: " << zLevelIndex << ", face " << Mesh::getString(Mesh::face_t(face)) << ", startPoint: " << startPoint << ", endPoint: " << endPoint << ", currentZ: " << currentZ << ", nBorderPointsXNew_: " << nBorderPointsXNew_;

        // determine boundary points between startPoint and endPoint, with same zLevel value
        // call stl_create_rings.create_ring_section_mesh

        PyObject *startPointPy = PythonUtility::convertToPython<Vec3>::get(startPoint);
        PyObject *endPointPy = PythonUtility::convertToPython<Vec3>::get(endPoint);
        //PyObject *loopSectionPy = PyObject_CallFunction(functionCreateRingSectionMesh_, "O O O f i", stlMeshPy, startPointPy, endPointPy, currentZ, nBorderPointsXNew_);
        PyObject *loopSectionPy = PyObject_CallFunction(functionCreateRingSection_, "s O O f i", stlFilename_.c_str(), startPointPy, endPointPy, currentZ, nBorderPointsXNew_);
        PythonUtility::checkForError();
        //  create_ring_section(input_filename, start_point, end_point, z_value, n_points)
        assert(loopSectionPy);

        std::vector<Vec3> loopSection = PythonUtility::convertFromPython<std::vector<Vec3>>::get(loopSectionPy);
        if (loopSection.size() != nBorderPointsXNew_)
          LOG(ERROR) << "Ring section from python script contains only " << loopSection.size() << " points, " << nBorderPointsXNew_ << " requested.";
        assert(loopSection.size() == nBorderPointsXNew_);

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
        // x = points that are not set by streamlines, but need to be set by the border

        // determine affected subdomains
        int subdomainIndex0 = 0;
        int subdomainIndex1 = 0;
        int facePerpendicular0 = 0;
        int facePerpendicular1 = 0;
        int pointIndexEdge = 0;
        int pointIndexCorner = 0;

        if (face == (int)Mesh::face_t::face0Minus)
        {
          subdomainIndex0 = 0;
          subdomainIndex1 = 2;
          facePerpendicular0 = (int)Mesh::face_t::face1Plus;
          facePerpendicular1 = (int)Mesh::face_t::face1Minus;
          pointIndexEdge = 0;
          pointIndexCorner = 0;
        }
        else if (face == (int)Mesh::face_t::face0Plus)
        {
          subdomainIndex0 = 1;
          subdomainIndex1 = 3;
          facePerpendicular0 = (int)Mesh::face_t::face1Plus;
          facePerpendicular1 = (int)Mesh::face_t::face1Minus;
          pointIndexEdge = nBorderPointsX_-1;
          pointIndexCorner = nBorderPointsXNew_-1;
        }
        else if (face == (int)Mesh::face_t::face1Minus)
        {
          subdomainIndex0 = 0;
          subdomainIndex1 = 1;
          facePerpendicular0 = (int)Mesh::face_t::face0Plus;
          facePerpendicular1 = (int)Mesh::face_t::face0Minus;
          pointIndexEdge = 0;
          pointIndexCorner = 0;
        }
        else if (face == (int)Mesh::face_t::face1Plus)
        {
          subdomainIndex0 = 2;
          subdomainIndex1 = 3;
          facePerpendicular0 = (int)Mesh::face_t::face0Plus;
          facePerpendicular1 = (int)Mesh::face_t::face0Minus;
          pointIndexEdge = nBorderPointsX_-1;
          pointIndexCorner = nBorderPointsXNew_-1;
        }

        // bottom subdomain
        if (zLevelIndex <= nBorderPointsZ_-1)
        {
          int zLevelIndexSubdomain = zLevelIndex;
          assert(zLevelIndexSubdomain < nBorderPointsZ_);
          //LOG(DEBUG) << "save for subdomains " << subdomainIndex0 << " and " << subdomainIndex1 << ", zLevelIndex " << zLevelIndexSubdomain;
          int i = 0;
          for (std::vector<Vec3>::const_iterator loopPointIter = loopSection.begin(); loopPointIter != loopSection.begin()+nBorderPointsX_; loopPointIter++, i++)
          {
            /*if (face == 1 && subdomainIndex0==1 && i == 0)
            {
              LOG(DEBUG) << "set borderPointsSubdomain[" << subdomainIndex0 << "][" << face << "][" << zLevelIndexSubdomain << "][" << i << "] to " << *loopPointIter;
            }*/

            assert(i < nBorderPointsX_);
            borderPointsSubdomain[subdomainIndex0][face][zLevelIndexSubdomain][i] = *loopPointIter;
          }

          i = 0;
          for (std::vector<Vec3>::const_iterator loopPointIter = loopSection.begin()+nBorderPointsX_-1; loopPointIter != loopSection.end(); loopPointIter++, i++)
          {
            assert(i < nBorderPointsX_);
            borderPointsSubdomain[subdomainIndex1][face][zLevelIndexSubdomain][i] = *loopPointIter;
          }

          // set end point of line of inner cross that does not have a traced streamline there, because it is at the border ("x" in figure)
          borderPointsSubdomain[subdomainIndex0][facePerpendicular0][zLevelIndexSubdomain][pointIndexEdge] = loopSection[nBorderPointsX_-1];
          borderPointsSubdomain[subdomainIndex1][facePerpendicular1][zLevelIndexSubdomain][pointIndexEdge] = loopSection[nBorderPointsX_-1];

          // set neighbouring border point at the end ("*" in figure)
          borderPointsSubdomain[subdomainIndex0][facePerpendicular1][zLevelIndexSubdomain][pointIndexCorner] = loopSection[0];
          borderPointsSubdomain[subdomainIndex1][facePerpendicular0][zLevelIndexSubdomain][pointIndexCorner] = loopSection[nBorderPointsXNew_-1];

#ifndef NDEBUG
#if 0
#ifdef STL_OUTPUT
          s.str("");
          s << "08_border_points_filled_face_" << Mesh::getString((Mesh::face_t)face) << "_z_" << zLevelIndexSubdomain << "_subdomain_" << subdomainIndex0;
          PyObject_CallFunction(functionOutputPoints_, "s i O f", s.str().c_str(), currentRankSubset_->ownRankNo(),
                                PythonUtility::convertToPython<std::vector<Vec3>>::get(borderPointsSubdomain[subdomainIndex0][face][zLevelIndexSubdomain]), 0.2);
          PythonUtility::checkForError();
          s.str("");
          s << "08_border_points_filled_face_" << Mesh::getString((Mesh::face_t)face) << "_z_" << zLevelIndexSubdomain << "_subdomain_" << subdomainIndex1;
          PyObject_CallFunction(functionOutputPoints_, "s i O f", s.str().c_str(), currentRankSubset_->ownRankNo(),
                                PythonUtility::convertToPython<std::vector<Vec3>>::get(borderPointsSubdomain[subdomainIndex1][face][zLevelIndexSubdomain]), 0.2);
          PythonUtility::checkForError();
#endif
#endif
#endif
          //std::copy(loopSection.begin(), loopSection.begin()+nBorderPointsX_, borderPointsSubdomain[subdomainIndex0][face][zLevelIndexSubdomain].begin());
          //std::copy(loopSection.begin()+nBorderPointsX_-1, loopSection.end(), borderPointsSubdomain[subdomainIndex1][face][zLevelIndexSubdomain].begin());
        }

        // top subdomain
        if (zLevelIndex >= nBorderPointsZ_-1)
        {
          subdomainIndex0 += 4;
          subdomainIndex1 += 4;
          int zLevelIndexSubdomain = zLevelIndex - (nBorderPointsZ_-1);
          assert(zLevelIndexSubdomain < nBorderPointsZ_);

          //LOG(DEBUG) << "save for subdomains " << subdomainIndex0 << " and " << subdomainIndex1 << ", zLevelIndex " << zLevelIndexSubdomain;
          //std::copy(loopSection.begin(), loopSection.begin()+nBorderPointsX_, borderPointsSubdomain[subdomainIndex0][face][zLevelIndexSubdomain].begin());
          //std::copy(loopSection.begin()+nBorderPointsX_-1, loopSection.end(), borderPointsSubdomain[subdomainIndex1][face][zLevelIndexSubdomain].begin());

          int i = 0;
          for (std::vector<Vec3>::const_iterator loopPointIter = loopSection.begin(); loopPointIter != loopSection.begin()+nBorderPointsX_; loopPointIter++, i++)
          {
            assert(i < nBorderPointsX_);
            borderPointsSubdomain[subdomainIndex0][face][zLevelIndexSubdomain][i] = *loopPointIter;
          }

          i = 0;
          for (std::vector<Vec3>::const_iterator loopPointIter = loopSection.begin()+nBorderPointsX_-1; loopPointIter != loopSection.end(); loopPointIter++, i++)
          {
            assert(i < nBorderPointsX_);
            borderPointsSubdomain[subdomainIndex1][face][zLevelIndexSubdomain][i] = *loopPointIter;
          }

          // set "center point" of loop section ("x" in figure)
          borderPointsSubdomain[subdomainIndex0][facePerpendicular0][zLevelIndexSubdomain][pointIndexEdge] = loopSection[nBorderPointsX_-1];
          borderPointsSubdomain[subdomainIndex1][facePerpendicular1][zLevelIndexSubdomain][pointIndexEdge] = loopSection[nBorderPointsX_-1];

          // set neighbouring border point at the end ("*" in figure)
          borderPointsSubdomain[subdomainIndex0][facePerpendicular1][zLevelIndexSubdomain][pointIndexCorner] = loopSection[0];
          borderPointsSubdomain[subdomainIndex1][facePerpendicular0][zLevelIndexSubdomain][pointIndexCorner] = loopSection[nBorderPointsXNew_-1];

#ifndef NDEBUG
#if 0
#ifdef STL_OUTPUT
          s.str("");
          s << "08_border_points_filled_face_" << Mesh::getString((Mesh::face_t)face) << "_z_" << zLevelIndexSubdomain << "_subdomain_" << subdomainIndex0;
          PyObject_CallFunction(functionOutputPoints_, "s i O f", s.str().c_str(), currentRankSubset_->ownRankNo(),
                                PythonUtility::convertToPython<std::vector<Vec3>>::get(borderPointsSubdomain[subdomainIndex0][face][zLevelIndexSubdomain]), 0.2);
          PythonUtility::checkForError();
          s.str("");
          s << "08_border_points_filled_face_" << Mesh::getString((Mesh::face_t)face) << "_z_" << zLevelIndexSubdomain << "_subdomain_" << subdomainIndex1;
          PyObject_CallFunction(functionOutputPoints_, "s i O f", s.str().c_str(), currentRankSubset_->ownRankNo(),
                                PythonUtility::convertToPython<std::vector<Vec3>>::get(borderPointsSubdomain[subdomainIndex1][face][zLevelIndexSubdomain]), 0.2);
          PythonUtility::checkForError();
#endif
#endif
#endif
        }

      }
    }  // for
  }

#ifndef NDEBUG
#ifdef STL_OUTPUT
#ifdef STL_OUTPUT_VERBOSE
  for (int subdomainIndex = 0; subdomainIndex < 8; subdomainIndex++)
  {
    std::stringstream s;
    s << "08_border_points_filled_subdomain_" << subdomainIndex;
    PyObject_CallFunction(functionOutputBorderPoints_, "s i O f", s.str().c_str(), currentRankSubset_->ownRankNo(),
                          PythonUtility::convertToPython<std::array<std::vector<std::vector<Vec3>>,4>>::get(borderPointsSubdomain[subdomainIndex]), 0.1);
    PythonUtility::checkForError();
  }
#endif
#endif
#endif

  // output points for debugging
#ifndef NDEBUG
  std::ofstream file("points.csv", std::ios::out | std::ios::trunc);
  assert (file.is_open());
  file << "subdomainIndex;faceNo;zLevelIndex;;points" << std::endl;
  for (int subdomainIndex = 0; subdomainIndex < 8; subdomainIndex++)
  {
    for (int faceNo = (int)Mesh::face_t::face0Minus; faceNo <= (int)Mesh::face_t::face1Plus; faceNo++)
    {
      for (int zLevelIndex = 0; zLevelIndex < borderPointsSubdomain[subdomainIndex][faceNo].size(); zLevelIndex++)
      {
        file << subdomainIndex << ";" << faceNo << ";" << zLevelIndex << ";;";
        for (int pointIndex = 0; pointIndex < borderPointsSubdomain[subdomainIndex][faceNo][zLevelIndex].size(); pointIndex++)
        {
          file << borderPointsSubdomain[subdomainIndex][faceNo][zLevelIndex][pointIndex][0] << ";" << borderPointsSubdomain[subdomainIndex][faceNo][zLevelIndex][pointIndex][1] << ";" << borderPointsSubdomain[subdomainIndex][faceNo][zLevelIndex][pointIndex][2] << ";";
        }
        file << std::endl;
      }
    }
  }
  file.close();
#endif
}

};  // namespace
