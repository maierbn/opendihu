#pragma once

#include <Python.h>  // has to be the first included header
#include <vector>
#include <streambuf>

#include <petscmat.h>
#include "control/types.h"

#include "utility/matrix.h"
#include "utility/matrix_operators.h"
#include "utility/vector_operators.h"
#include "utility/vector_operators_multiplication.h"

namespace MathUtility
{

//! returns v*v (square)
inline double sqr(double v){return v*v;}

//! returns v*v (square)
inline Vc::double_v sqr(Vc::double_v v){return v*v;}

//! returns v*v (square)
inline int sqr(int v){return v*v;}

template<typename double_v_t>
double_v_t pow(double_v_t base, double exponent);

//! return the euclidean norm of the difference vector, i.e. the distance between the nodes
template<int D>
double distance(const VecD<D> node1, const VecD<D> node2);

//! return the euclidean norm of the vector, 1D vector
template<typename double_v_t=double>
double_v_t length(const VecD<1,double_v_t> node);

//! return the euclidean norm of the vector, 2D vector
template<typename double_v_t=double>
double_v_t length(const VecD<2,double_v_t> node);

//! return the euclidean norm of the vector, 3D vector
template<typename double_v_t=double>
double_v_t length(const VecD<3,double_v_t> node);

//! same as length, return the euclidean norm of the vector
template<int D, typename double_v_t=double>
double_v_t norm(const VecD<D,double_v_t> node);

//! compute the squared norm of the vector
template<int D>
double normSquared(const VecD<D> node);

//! compute the squared norm of the vector
template<int D>
double normSquared(const VecD<D,Vc::double_v> node);

//! return the normalized vector
template<int D, typename double_v_t=double>
VecD<D> normalized(VecD<D,double_v_t> &vector);

//! normalize the vector
template<int D, typename double_v_t=double>
void normalize(VecD<D,double_v_t> &vector);

//! arc cosine
template<typename double_v_t=double>
double_v_t acos(double_v_t value);

//! fabs for double or Vc::double_v
template<typename double_v_t=double>
double_v_t abs(double_v_t value);

//! compute a value h that depends linearly on the mesh width, needed as factor for consistent regularization (i.e. ε→0 for h→0)
template<typename double_v_t, int nNodes>
double_v_t computeApproximateMeshWidth(const std::array<VecD<3,double_v_t>,nNodes> &geometryValues);

//! compute the matrix T = J^{-1}J^{-T} when J is the jacobian
//! returns the matrix in row major storage order, matrix is symmetric and the determinant of the jacobian
template<typename double_v_t=double>
std::array<double_v_t,9> computeTransformationMatrixAndDeterminant(const std::array<VecD<3,double_v_t>,3> &jacobian, double_v_t &determinant);

//! compute the matrix T = J^{-1} A J^{-T} when J is the jacobian and A the symmetric diffusionTensor
//! returns the matrix in row major storage order, matrix is symmetric and the determinant of the jacobian
template<typename double_v_t=double>
std::array<double_v_t,9> computeTransformationDiffusionMatrixAndDeterminant(const std::array<VecD<3,double_v_t>,3> &jacobian, const Matrix<3,3,double_v_t> &diffusionTensor, double_v_t &determinant);

//! computes the determinant of the matrix, for 2D matrix
template<typename double_v_t=double>
double_v_t computeDeterminant(const Tensor2<2,double_v_t> &jacobian);

//! computes the determinant of the matrix, for 3D matrix
template<typename double_v_t=double>
double_v_t computeDeterminant(const Tensor2<3,double_v_t> &jacobian);

//! computes the inverse of the symmetric matrix and the determinant, for 2D matrix
template<typename double_v_t=double>
Tensor2<2,double_v_t> computeSymmetricInverse(const Tensor2<2,double_v_t> &matrix, double_v_t &approximateMeshWidth, double_v_t &determinant);

//! computes the inverse of the symmetric matrix and the determinant, for 3D matrix
template<typename double_v_t=double>
Tensor2<3,double_v_t> computeSymmetricInverse(const Tensor2<3,double_v_t> &matrix, double_v_t &approximateMeshWidth, double_v_t &determinant);

//! computes the inverse of the non-symmetric matrix and the determinant, for 1D matrix
template<typename double_v_t=double>
Tensor2<1,double_v_t> computeInverse(const Tensor2<1,double_v_t> &matrix, double_v_t &approximateMeshWidth, double_v_t &determinant);

//! computes the inverse of the non-symmetric matrix and the determinant, for 2D matrix
template<typename double_v_t=double>
Tensor2<2,double_v_t> computeInverse(const Tensor2<2,double_v_t> &matrix, double_v_t &approximateMeshWidth, double_v_t &determinant);

//! computes the inverse of the non-symmetric matrix and the determinant, for 3D matrix
template<typename double_v_t=double>
Tensor2<3,double_v_t> computeInverse(const Tensor2<3,double_v_t> &matrix, double_v_t &approximateMeshWidth, double_v_t &determinant);

//! computes the cofactor matrix of the non-symmetric matrix. inv = 1/det * adj, adj = cof^T, i.e. cof(M) = det(M) * M^{-T}
template<typename double_v_t>
Tensor2<2,double_v_t> computeCofactorMatrix(const Tensor2<2,double_v_t> &matrix);

template<typename double_v_t>
Tensor2<3,double_v_t> computeCofactorMatrix(const Tensor2<3,double_v_t> &matrix);

//! computes the transpose matrix
template<typename double_v_t>
Tensor2<3,double_v_t> computeTranspose(const Tensor2<3,double_v_t> &matrix);

//! transform a 3xD2 matrix to a DxD matrix by filling up with identity entries
template<int D, int D2>
std::array<std::array<double,D>,D> transformToDxD(const std::array<Vec3,D2> &matrix);

//! transform a 3xD2 matrix to a DxD matrix by filling up with identity entries, vectorized version
template<int D, int D2>
std::array<std::array<Vc::double_v,D>,D> transformToDxD(const std::array<VecD<3,Vc::double_v>,D2> &matrix);

//! transform a vector with D2 entries to a vector with D entries
template<int D, int D2>
VecD<D> transformToD(const VecD<D2> &vector);

//! computes v1^T * T * v2 where T is the symmetric transformation matrix (3D)
template<typename double_v1_t=double, typename double_v2_t=double, typename double_v3_t=double>
double_v1_t applyTransformation(const std::array<double_v1_t,9> &transformationMatrix, const VecD<3,double_v2_t> &vector1, const VecD<3,double_v3_t> &vector2);

//! computes v1^T * T * v2 where T is the symmetric transformation matrix (2D)
template<typename double_v1_t=double, typename double_v2_t=double, typename double_v3_t=double>
double_v1_t applyTransformation(const std::array<double_v1_t,4> &transformationMatrix, const VecD<2,double_v2_t> &vector1, const VecD<2,double_v3_t> &vector2);

//! computes the factor J_D that is needed when transforming an integral from world coordinates to parameter space
template<typename double_v_t=double>
double_v_t computeIntegrationFactor(const std::array<VecD<3,double_v_t>,1> &jacobian);

template<typename double_v_t=double>
double_v_t computeIntegrationFactor(const std::array<VecD<3,double_v_t>,2> &jacobian);

template<typename double_v_t=double>
double_v_t computeIntegrationFactor(const std::array<VecD<3,double_v_t>,3> &jacobian);

//! check if the value is not nan or inf
template<typename double_v_t=double>
bool isFinite(double_v_t value);

//! compute 3D cross product
Vec3 cross(const Vec3 &vector1, const Vec3 &vector2);

//! compute dot product
double dot(const Vec3 &vector1, const Vec3 &vector2);

//! return a values of the Levi-Civita permutation symbol
int permutation(int i, int j, int k);

const double PI = 3.14159265358979323846;

//! check if vector b is a subsequence of a, i.e. is contained in a
bool isSubsequenceOf(std::vector<int> a, std::vector<int> b, size_t &subsequenceAStartPos);

//! check if the two vectors have equal entries
template<int D>
bool equals(std::array<double,D> a, std::array<double,D> b, double tolerance=1e-15);

//! rotate the matrix such that unit vector (1,0,0) now points to directionVector, directionVector does not need to be normalized
template<typename double_v_t>
void rotateMatrix(Matrix<2,2,double_v_t> &matrix, VecD<2,double_v_t> directionVector);
//! rotate the matrix such that unit vector (1,0,0) now points to directionVector, directionVector does not need to be normalized

template<typename double_v_t>
void rotateMatrix(Matrix<3,3,double_v_t> &matrix, VecD<3,double_v_t> directionVector);

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

//! check if the given value contains any value > 1e+75 or nan
bool containsNanOrInf(const double value);

//! check if the given value contains any value > 1e+75 or nan
bool containsNanOrInf(const Vc::double_v value);

template<typename T, std::size_t N>
bool containsNanOrInf(std::array<T,N> value);

/** pow as constexpr function
 *  https://stackoverflow.com/a/27270730/10290071
 */
template<typename T, typename U>
T constexpr powConst(T base, U exponent)
{
  static_assert(std::is_integral<U>(), "exponent must be integral");
  return exponent == 0 ? 1 : base * powConst(base, exponent - 1);
}

}  // namespace

#include "utility/math_utility.tpp"
