#pragma once

#include <Python.h>  // has to be the first included header
#include <vector>
#include <streambuf>

#include <petscmat.h>
#include "control/types.h"
#include "semt/Semt.h"    // include semt for operator<<(std::vector<double>)

#include "utility/matrix.h"
#include "utility/matrix_operators.h"
#include "utility/vector_operators.h"

namespace MathUtility 
{

//! returns v*v (square)
double sqr(double v);

//! returns v*v (square)
int sqr(int v);

//! return the euclidean norm of the difference vector, i.e. the distance between the nodes
double distance(const Vec3 node1, const Vec3 node2);

//! return the euclidean norm of the vector
double length(const Vec3 node);

//! same as length, return the euclidean norm of the vector
double norm(const Vec3 node);

//! compute the matrix T = J^{-1}J^{-T} when J is the jacobian
//! returns the matrix in row major storage order, matrix is symmetric and the determinant of the jacobian
std::array<double,9> computeTransformationMatrixAndDeterminant(const std::array<Vec3,3> &jacobian, double &determinant);

//! computes the determinant of the matrix 
double computeDeterminant(const std::array<Vec3,3> &jacobian);

//! computes the inverse of the symmetric matrix and the determinant
std::array<Vec3,3> computeSymmetricInverse(const std::array<Vec3,3> &matrix, double &determinant);

//! computes the inverse of the non-symmetric matrix and the determinant
std::array<Vec3,3> computeInverse(const std::array<Vec3,3> &matrix, double &determinant);

//! computes v1^T * T * v2 where T is the symmetric transformation matrix (3D)
double applyTransformation(const std::array<double,9> &transformationMatrix, const Vec3 &vector1, const Vec3 &vector2);

//! computes v1^T * T * v2 where T is the symmetric transformation matrix (2D)
double applyTransformation(const std::array<double,4> &transformationMatrix, const Vec2 &vector1, const Vec2 &vector2);  

//! computes the factor J_D that is needed when transforming an integral from world coordinates to parameter space
template<int D>
double computeIntegrationFactor(const std::array<Vec3,D> &jacobian);

//! return a values of the Levi-Civita permutation symbol
int permutation(int i, int j, int k);

const double PI = 3.14159265358979323846;

//! check if vector b is a subsequence of a, i.e. is contained in a
bool isSubsequenceOf(std::vector<int> a, std::vector<int> b, size_t &subsequenceAStartPos);
 
}  // namespace
