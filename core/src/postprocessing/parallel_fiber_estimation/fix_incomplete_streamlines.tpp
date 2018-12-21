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

  LOG(DEBUG) << "valid streamlines:";
  // output which streamlines are valid
  for (int subdomainIndex = 0; subdomainIndex < 8; subdomainIndex++)
  {
    for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
    {
      for (int pointIndex = 0; pointIndex < nBorderPointsX_; pointIndex++)
      {
        LOG(DEBUG) << "subdomain " << subdomainIndex << ", face " << Mesh::getString((Mesh::face_t)face) << ", streamline " << pointIndex << ": "
          << std::boolalpha << borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Minus][pointIndex];
      }
    }
  }

  // fill invalid streamlines
  for (int subdomainIndex = 0; subdomainIndex < 8; subdomainIndex++)
  {
    for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
    {
      int lastValid = 0;
      bool lastStreamlinesWereInvalid = false;
      for (int pointIndex = 0; pointIndex < nBorderPointsX_; pointIndex++)
      {
        // if current streamline is valid
        if (borderPointsSubdomainAreValid[subdomainIndex][face][pointIndex])
        {
          if (lastStreamlinesWereInvalid)
          {
            // loop over streamlines between lastValid and pointIndex, these are all invalid, interpolate them
            for (int invalidStreamlineIndex = lastValid+1; invalidStreamlineIndex < pointIndex; invalidStreamlineIndex++)
            {
              int zLevelIndex = 0;
              if (streamlineDirectionUpwards)
              {
                zLevelIndex = 0;
              }
              else
              {
                zLevelIndex = nBorderPointsZNew_-1;
              }

              double alpha = MathUtility::norm<3>(borderPointsSubdomain[subdomainIndex][face][0][invalidStreamlineIndex] - borderPointsSubdomain[subdomainIndex][face][zLevelIndex][lastValid])
                / MathUtility::norm<3>(borderPointsSubdomain[subdomainIndex][face][zLevelIndex][pointIndex] - borderPointsSubdomain[subdomainIndex][face][zLevelIndex][lastValid]);

              //loop over points of streamline from bottom to top
              for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZNew_; zLevelIndex++)
              {
                borderPointsSubdomain[subdomainIndex][face][zLevelIndex][invalidStreamlineIndex]
                  = (1.-alpha) * borderPointsSubdomain[subdomainIndex][face][zLevelIndex][lastValid]
                    + alpha * borderPointsSubdomain[subdomainIndex][face][zLevelIndex][pointIndex];
              }
              LOG(DEBUG) << "interpolate in subdomain " << subdomainIndex << ", face " << Mesh::getString((Mesh::face_t)face)
                << " streamline " << invalidStreamlineIndex << " from " << lastValid << " and " << pointIndex << ", alpha: " << alpha;
              LOG(DEBUG) << "invalid seed: " << borderPointsSubdomain[subdomainIndex][face][0][invalidStreamlineIndex] << ", lastValid start: " << borderPointsSubdomain[subdomainIndex][face][zLevelIndex][lastValid]
                << " next valid: " << borderPointsSubdomain[subdomainIndex][face][zLevelIndex][pointIndex];
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
