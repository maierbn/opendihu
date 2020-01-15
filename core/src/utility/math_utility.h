#pragma once

#include <Python.h>  // has to be the first included header
#include <vector>
#include <streambuf>

#include <petscmat.h>
#include "control/types.h"

#include "utility/matrix.h"
#include "utility/matrix_operators.h"
#include "utility/vector_operators.h"

namespace MathUtility
{

//! returns v*v (square)
inline double sqr(double v){return v*v;}

//! returns v*v (square)
inline int sqr(int v){return v*v;}

//! return the euclidean norm of the difference vector, i.e. the distance between the nodes
template<int D>
double distance(const VecD<D> node1, const VecD<D> node2);

//! return the euclidean norm of the vector
template<int D>
double length(const VecD<D> node);

//! same as length, return the euclidean norm of the vector
template<int D>
double norm(const VecD<D> node);

//! compute the squared norm of the vector
template<int D>
double normSquared(const VecD<D> node);

//! return the normalized vector
template<int D>
VecD<D> normalized(VecD<D> &vector);

//! normalize the vector
template<int D>
void normalize(VecD<D> &vector);

//! compute the matrix T = J^{-1}J^{-T} when J is the jacobian
//! returns the matrix in row major storage order, matrix is symmetric and the determinant of the jacobian
std::array<double,9> computeTransformationMatrixAndDeterminant(const std::array<Vec3,3> &jacobian, double &determinant);

//! compute the matrix T = J^{-1} A J^{-T} when J is the jacobian and A the symmetric diffusionTensor
//! returns the matrix in row major storage order, matrix is symmetric and the determinant of the jacobian
std::array<double,9> computeTransformationDiffusionMatrixAndDeterminant(const std::array<Vec3,3> &jacobian, const Matrix<3,3> &diffusionTensor, double &determinant);

//! computes the determinant of the matrix
template<int D>
double computeDeterminant(const Tensor2<D> &jacobian);

//! computes the inverse of the symmetric matrix and the determinant
template<int D>
Tensor2<D> computeSymmetricInverse(const Tensor2<D> &matrix, double &determinant);

//! computes the inverse of the non-symmetric matrix and the determinant
template<int D>
Tensor2<D> computeInverse(const Tensor2<D> &matrix, double &determinant);

//! computes the cofactor matrix of the non-symmetric matrix. inv = 1/det * adj, adj = cof^T, i.e. cof(M) = det(M) * M^{-T}
template<int D>
Tensor2<D> computeCofactorMatrix(const Tensor2<D> &matrix);

//! transform a 3xD2 matrix to a DxD matrix by filling up with identity entries
template<int D, int D2>
std::array<std::array<double,D>,D> transformToDxD(const std::array<Vec3,D2> &matrix);

//! transform a vector with D2 entries to a vector with D entries
template<int D, int D2>
VecD<D> transformToD(const VecD<D2> &vector);

//! computes v1^T * T * v2 where T is the symmetric transformation matrix (3D)
double applyTransformation(const std::array<double,9> &transformationMatrix, const Vec3 &vector1, const Vec3 &vector2);

//! computes v1^T * T * v2 where T is the symmetric transformation matrix (2D)
double applyTransformation(const std::array<double,4> &transformationMatrix, const Vec2 &vector1, const Vec2 &vector2);

//! computes the factor J_D that is needed when transforming an integral from world coordinates to parameter space
template<int D>
double computeIntegrationFactor(const std::array<Vec3,D> &jacobian);

//! compute 3D cross product
Vec3 cross(const Vec3 &vector1, const Vec3 &vector2);

//! compute dot product
double dot(const Vec3 &vector1, const Vec3 &vector2);

//! return a values of the Levi-Civita permutation symbol
int permutation(int i, int j, int k);

const double PI = 3.14159265358979323846;

//! check if vector b is a subsequence of a, i.e. is contained in a
bool isSubsequenceOf(std::vector<int> a, std::vector<int> b, size_t &subsequenceAStartPos);

//! rotate the matrix such that unit vector (1,0,0) now points to directionVector, directionVector does not need to be normalized
template<int D>
void rotateMatrix(Matrix<D,D> &matrix, Vec3 directionVector);

//! compute the inverse Phi mapping, i.e. get the xi0, xi1 coordinates of a point in a quadrilateral
void quadrilateralGetPointCoordinates(const std::array<Vec3,4> quadrilateral, const Vec3 point, Vec2 &xi);

//! read a Vec3 from an open fstream at the current position of the file get pointer
template<typename T>
void readPoint(T &file, Vec3 &point);

//! write a Vec3 to an open fstream at the current position of the file put pointer
template<typename T>
void writePoint(T &file, Vec3 &point);

//! estimate the maximum eigenvalue of the matrix by using the power iteration / von Mises iteration algorithm, the algorithm is stopped after 15 iterations
double estimateMaximumEigenvalue(const Tensor2<3> &matrix);

//! estimate the condition number |sigma_max/sigma_min| of the matrix
double estimateConditionNumber(const Tensor2<3> &matrix, const Tensor2<3> &inverseMatrix);

/** pow as constexpr function
 *  https://stackoverflow.com/a/27270730/10290071
 */
template<typename T, typename U>
T constexpr pow(T base, U exponent)
{
  static_assert(std::is_integral<U>(), "exponent must be integral");
  return exponent == 0 ? 1 : base * pow(base, exponent - 1);
}

}  // namespace

#include "utility/math_utility.tpp"
