#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
refineBoundaryPoints(std::array<std::vector<std::vector<Vec3>>,4> &boundaryPointsOld, std::array<std::vector<std::vector<Vec3>>,4> &boundaryPoints)
{
  // refine boundaryPoints to twice the precision, only in x and y direction, stays the same in z direction
  // then we have nBoundaryPointsXNew_ points per x,y-direction and nBoundaryPointsZ_ in z direction, each time also including first and last boundary point
  // boundaryPoints[face_t][z-level][pointIndex]
  for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
  {

    boundaryPoints[face].resize(nBoundaryPointsZ_);
    for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZ_; zLevelIndex++)
    {
      boundaryPoints[face][zLevelIndex].resize(nBoundaryPointsXNew_);
      Vec3 previousPoint = boundaryPointsOld[face][zLevelIndex][0];

      boundaryPoints[face][zLevelIndex][0] = previousPoint;

      for (int pointIndex = 1; pointIndex < nBoundaryPointsX_; pointIndex++)
      {
        const Vec3 &currentPoint = boundaryPointsOld[face][zLevelIndex][pointIndex];
        boundaryPoints[face][zLevelIndex][2*pointIndex-1] = 0.5*(currentPoint + previousPoint);

        // interpolate point in between
        boundaryPoints[face][zLevelIndex][2*pointIndex+0] = currentPoint;

        previousPoint = currentPoint;
      }
    }
  }

#ifndef NDEBUG
#ifdef STL_OUTPUT
#ifdef STL_OUTPUT_VERBOSE
  PyObject_CallFunction(functionOutputBoundaryPoints_, "s i i O f", "01_boundary_points_old", currentRankSubset_->ownRankNo(), level_,
                        PythonUtility::convertToPython<std::array<std::vector<std::vector<Vec3>>,4>>::get(boundaryPointsOld), 0.02);
  PythonUtility::checkForError();

  for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
  {
    std::stringstream s;
    s << "02_boundary_points_face_" << Mesh::getString((Mesh::face_t)face);
    PyObject_CallFunction(functionOutputPoints_, "s i i O f", s.str().c_str(), currentRankSubset_->ownRankNo(), level_,
                          PythonUtility::convertToPython<std::vector<Vec3>>::get(boundaryPoints[face][0]), 0.2);
    PythonUtility::checkForError();
  }

  PyObject_CallFunction(functionOutputBoundaryPoints_, "s i i O f", "02_boundary_points", currentRankSubset_->ownRankNo(), level_,
                        PythonUtility::convertToPython<std::array<std::vector<std::vector<Vec3>>,4>>::get(boundaryPoints), 0.03);
  PythonUtility::checkForError();
#endif
#endif
#endif
}

} // namespace
