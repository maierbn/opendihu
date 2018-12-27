#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
fixIncompleteStreamlines(std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &borderPointsSubdomain,
                         std::array<std::array<std::vector<bool>,4>,8> &borderPointsSubdomainAreValid, bool streamlineDirectionUpwards)
{
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
}

};  // namespace
