#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
fixIncompleteStreamlines(std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &borderPointsSubdomain,
                         std::array<std::array<std::vector<bool>,4>,8> &borderPointsSubdomainAreValid)
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

}

};  // namespace
