#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
refineBorderPoints(std::array<std::vector<std::vector<Vec3>>,4> &borderPointsOld, std::array<std::vector<std::vector<Vec3>>,4> &borderPoints)
{
  // refine borderPoints to twice the precision, only in x and y direction, stays the same in z direction
  // then we have nBorderPointsXNew_ points per x,y-direction and nBorderPointsZ_ in z direction, each time also including first and last border point
  // borderPoints[face_t][z-level][pointIndex]
  for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
  {

    borderPoints[face].resize(nBorderPointsZ_);
    for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
    {
      borderPoints[face][zLevelIndex].resize(nBorderPointsXNew_);
      Vec3 previousPoint = borderPointsOld[face][zLevelIndex][0];

      borderPoints[face][zLevelIndex][0] = previousPoint;

      for (int pointIndex = 1; pointIndex < nBorderPointsX_; pointIndex++)
      {
        const Vec3 &currentPoint = borderPointsOld[face][zLevelIndex][pointIndex];
        borderPoints[face][zLevelIndex][2*pointIndex-1] = 0.5*(currentPoint + previousPoint);

        // interpolate point in between
        borderPoints[face][zLevelIndex][2*pointIndex+0] = currentPoint;

        previousPoint = currentPoint;
      }
    }
  }

#ifndef NDEBUG
  PyObject_CallFunction(functionOutputBorderPoints_, "s i O f", "01_border_points_old", currentRankSubset_->ownRankNo(),
                        PythonUtility::convertToPython<std::array<std::vector<std::vector<Vec3>>,4>>::get(borderPointsOld), 0.2);
  PythonUtility::checkForError();

  for (int face = Mesh::face_t::face0Minus; face <= Mesh::face_t::face1Plus; face++)
  {
    std::stringstream s;
    s << "02_border_points_face_" << Mesh::getString((Mesh::face_t)face);
    PyObject_CallFunction(functionOutputPoints_, "s i O f", s.str().c_str(), currentRankSubset_->ownRankNo(),
                          PythonUtility::convertToPython<std::vector<Vec3>>::get(borderPoints[face][0]), 0.2);
    PythonUtility::checkForError();
  }

  PyObject_CallFunction(functionOutputBorderPoints_, "s i O f", "02_border_points", currentRankSubset_->ownRankNo(),
                        PythonUtility::convertToPython<std::array<std::vector<std::vector<Vec3>>,4>>::get(borderPoints), 0.3);
  PythonUtility::checkForError();

#endif
}

};  // namespace
