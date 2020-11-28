#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
outputBoundaryPoints(std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &boundaryPointsSubdomain, std::string name)
{
#ifndef NDEBUG
  std::vector<Vec3> points;

  LOG(DEBUG) << "outputBoundaryPoints(name=" << name << ")";

  // boundaryPointsSubdomain[subdomain index][face_t][z-level][point index]

  for (int subdomainIndex = 0; subdomainIndex < 8; subdomainIndex++)
  {
    for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
    {
      for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZ_; zLevelIndex++)
      {
        for (int pointIndex = 0; pointIndex != boundaryPointsSubdomain[subdomainIndex][face][zLevelIndex].size(); pointIndex++)
        {
          points.push_back(boundaryPointsSubdomain[subdomainIndex][face][zLevelIndex][pointIndex]);
        }
      }
    }
  }

  std::stringstream filename;
  filename << name;
  PyObject_CallFunction(functionOutputPoints_, "s i i O f", filename.str().c_str(), currentRankSubset_->ownRankNo(), level_,
                        PythonUtility::convertToPython<std::vector<Vec3>>::get(points), 0.1);
  PythonUtility::checkForError();
#endif
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
outputStreamlines(std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &boundaryPointsSubdomain, std::string name)
{
#ifndef NDEBUG
  std::vector<std::vector<Vec3>> streamlines;

  // boundaryPointsSubdomain[subdomain index][face_t][z-level][point index]

  //LOG(DEBUG) << "outputStreamlines: ";
  for (int subdomainIndex = 0; subdomainIndex < 8; subdomainIndex++)
  {
    for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
    {
      //LOG(DEBUG) << "subdomainIndex: " << subdomainIndex << ", face: " << Mesh::getString((Mesh::face_t)face);

      for (int pointIndex = 0; pointIndex != boundaryPointsSubdomain[subdomainIndex][face][0].size(); pointIndex++)
      {
        std::vector<Vec3> points;
        for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZ_; zLevelIndex++)
        {
          points.push_back(boundaryPointsSubdomain[subdomainIndex][face][zLevelIndex][pointIndex]);
        }
        streamlines.push_back(points);
      }
    }
  }

  std::stringstream filename;
  filename << name;
  PyObject_CallFunction(functionOutputStreamlines_, "s i i O f", filename.str().c_str(), currentRankSubset_->ownRankNo(), level_,
                        PythonUtility::convertToPython<std::vector<std::vector<Vec3>>>::get(streamlines), 0.1);
  PythonUtility::checkForError();
#endif
}

} // namespace
