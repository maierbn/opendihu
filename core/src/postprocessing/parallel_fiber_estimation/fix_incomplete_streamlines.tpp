#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
fixIncompleteStreamlines(std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &borderPointsSubdomain,
                         std::array<std::array<std::vector<bool>,4>,8> &borderPointsSubdomainAreValid, bool streamlineDirectionUpwards,
                         std::array<std::vector<std::vector<Vec3>>,4> borderPoints)
{
  // std::array<std::vector<std::vector<Vec3>>,4> borderPoints;    // [face_t][z-level][pointIndex]
  // borderPointsSubdomain[subdomainIndex][face][zLevelIndex][pointIndex]
  // borderPointsSubdomainAreValid[subdomainIndex][face][pointIndex]

  //TODO fix streamlines at the corner of the domain frmo the border points if they ar e invalid

  goto fynn;
  fynn:

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


  //
  //   ^ --(1+)-> ^   ^ --(1+)-> ^
  //   0-   [2]   0+  0-   [3]   0+
  //   | --(1-)-> |   | --(1-)-> |
  //
  //   ^ --(1+)-> ^   ^ --(1+)-> ^
  // ^ 0-   [0]   0+  0-   [1]   0+
  // | | --(1-)-> |   | --(1-)-> |
  // +-->

  std::vector<std::tuple<int,int,int>> subdomainFacePointIndex {
    {0, Mesh::face_t::face0Minus, 0}, {0, Mesh::face_t::face1Minus, 0},
    {1, Mesh::face_t::face1Minus, nBorderPointsX_-1}, {1, Mesh::face_t::face0Plus, 0},
    {2, Mesh::face_t::face0Minus, nBorderPointsX_-1}, {1, Mesh::face_t::face1Plus, 0},
    {3, Mesh::face_t::face1Plus, nBorderPointsX_-1}, {1, Mesh::face_t::face0Plus, nBorderPointsX_-1}
  };

  // TODO: center point

  for (int i = 0; i < 8; i+=2)
  {
    int subdomainIndex0 = std::get<0>(subdomainFacePointIndex[i]);
    int face0 = std::get<1>(subdomainFacePointIndex[i]);
    int pointIndex0 = std::get<2>(subdomainFacePointIndex[i]);

    int subdomainIndex1 = std::get<0>(subdomainFacePointIndex[i+1]);
    int face1 = std::get<1>(subdomainFacePointIndex[i+1]);
    int pointIndex1 = std::get<2>(subdomainFacePointIndex[i+1]);

    // if at least one of the two instances of the corner streamline is invalid
    if (!borderPointsSubdomainAreValid[subdomainIndex0][face0][pointIndex0]
      || !borderPointsSubdomainAreValid[subdomainIndex1][face1][pointIndex1])
    {
      // derive the corner streamline from the initial border points

      LOG(DEBUG) << "subdomain " << subdomainIndex0 << " face " << Mesh::getString((Mesh::face_t)face0)
        << " pointIndex " << pointIndex0 << " valid? " << borderPointsSubdomainAreValid[subdomainIndex0][face0][pointIndex0];
      LOG(DEBUG) << "subdomain " << subdomainIndex1 << " face " << Mesh::getString((Mesh::face_t)face1)
        << " pointIndex " << pointIndex1 << " valid? " << borderPointsSubdomainAreValid[subdomainIndex1][face1][pointIndex1];

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

};  // namespace
