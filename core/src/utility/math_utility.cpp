#include "utility/math_utility.h"

#include <cmath>

namespace MathUtility
{
/*
double sqr(double v)
{
  return v*v;
}

int sqr(int v)
{
  return v*v;
}*/

template<>
double length<2>(const Vec2 node)
{
  return sqrt(sqr(node[0])
    + sqr(node[1]));
}

template<>
double length<3>(const Vec3 node)
{
  return sqrt(sqr(node[0])
    + sqr(node[1])
    + sqr(node[2]));
}

template<>
double distance<3>(const Vec3 node1, const Vec3 node2)
{
  return sqrt(sqr(node1[0]-node2[0])
    + sqr(node1[1]-node2[1])
    + sqr(node1[2]-node2[2]));
}

template<>
double distance<2>(const Vec2 node1, const Vec2 node2)
{
  return sqrt(sqr(node1[0]-node2[0])
    + sqr(node1[1]-node2[1]));
}

template<>
double normSquared<1>(const VecD<1> node)
{
  return sqr(node[0]);
}

template<>
double normSquared<2>(const VecD<2> node)
{
  return sqr(node[0]) + sqr(node[1]);
}

template<>
double normSquared<3>(const VecD<3> node)
{
  return sqr(node[0]) + sqr(node[1]) + sqr(node[2]);
}


std::array<double,9> computeTransformationMatrixAndDeterminant(const std::array<Vec3,3> &jacobian, double &determinant)
{
  // rename input values
  const double m11 = jacobian[0][0];
  const double m21 = jacobian[0][1];
  const double m31 = jacobian[0][2];
  const double m12 = jacobian[1][0];
  const double m22 = jacobian[1][1];
  const double m32 = jacobian[1][2];
  const double m13 = jacobian[2][0];
  const double m23 = jacobian[2][1];
  const double m33 = jacobian[2][2];

  // the following code is copy&pasted from the output of ./invert_mapping.py
  determinant = m11*m22*m33 - m11*m23*m32 - m12*m21*m33 + m12*m23*m31 + m13*m21*m32 - m13*m22*m31;
  double prefactor = 1./sqr(determinant);

  // 0 1 2
  // 3 4 5
  // 6 7 8

  // compute matrix entries of symmetric matrix
  std::array<double, 9> result;
  result[0] = prefactor * (sqr(m12*m23 - m13*m22) + sqr(m12*m33 - m13*m32) + sqr(m22*m33 - m23*m32));
  result[1] = prefactor * (-(m11*m23 - m13*m21)*(m12*m23 - m13*m22) - (m11*m33 - m13*m31)*(m12*m33 - m13*m32) - (m21*m33 - m23*m31)*(m22*m33 - m23*m32));
  result[2] = prefactor * ((m11*m22 - m12*m21)*(m12*m23 - m13*m22) + (m11*m32 - m12*m31)*(m12*m33 - m13*m32) + (m21*m32 - m22*m31)*(m22*m33 - m23*m32));
  result[3] = result[1];
  result[4] = prefactor * (sqr(m11*m23 - m13*m21) + sqr(m11*m33 - m13*m31) + sqr(m21*m33 - m23*m31));
  result[5] = prefactor * (-(m11*m22 - m12*m21)*(m11*m23 - m13*m21) - (m11*m32 - m12*m31)*(m11*m33 - m13*m31) - (m21*m32 - m22*m31)*(m21*m33 - m23*m31));
  result[6] = result[2];
  result[7] = result[5];
  result[8] = prefactor * (sqr(m11*m22 - m12*m21) + sqr(m11*m32 - m12*m31) + sqr(m21*m32 - m22*m31));

  return result;
}


std::array<double,9> computeTransformationDiffusionMatrixAndDeterminant(const std::array<Vec3,3> &jacobian, const Matrix<3,3> &diffusionTensor, double &determinant)
{
  // rename input values
  const double m11 = jacobian[0][0];
  const double m21 = jacobian[0][1];
  const double m31 = jacobian[0][2];
  const double m12 = jacobian[1][0];
  const double m22 = jacobian[1][1];
  const double m32 = jacobian[1][2];
  const double m13 = jacobian[2][0];
  const double m23 = jacobian[2][1];
  const double m33 = jacobian[2][2];

  const double a11 = diffusionTensor(0,0);
  const double a21 = diffusionTensor(0,1);
  const double a31 = diffusionTensor(0,2);
  const double a12 = diffusionTensor(1,0);
  const double a22 = diffusionTensor(1,1);
  const double a32 = diffusionTensor(1,2);
  const double a13 = diffusionTensor(2,0);
  const double a23 = diffusionTensor(2,1);
  const double a33 = diffusionTensor(2,2);

  // the following code is copy&pasted from the output of ./invert_mapping.py
  determinant = m11*m22*m33 - m11*m23*m32 - m12*m21*m33 + m12*m23*m31 + m13*m21*m32 - m13*m22*m31;
  double prefactor = 1./sqr(determinant);

  // 0 1 2
  // 3 4 5
  // 6 7 8

  // compute matrix entries of symmetric matrix  J^-1 A J^-T  assuming that A is symmetric
  std::array<double, 9> result;
  result[0] = prefactor * ((m12*m23 - m13*m22)*(a13*(m22*m33 - m23*m32) - a23*(m12*m33 - m13*m32) + a33*(m12*m23 - m13*m22)) - (m12*m33 - m13*m32)*(a12*(m22*m33 - m23*m32) - a22*(m12*m33 - m13*m32) + a32*(m12*m23 - m13*m22)) + (m22*m33 - m23*m32)*(a11*(m22*m33 - m23*m32) - a21*(m12*m33 - m13*m32) + a31*(m12*m23 - m13*m22)));
  result[1] = prefactor * (-(m11*m23 - m13*m21)*(a13*(m22*m33 - m23*m32) - a23*(m12*m33 - m13*m32) + a33*(m12*m23 - m13*m22)) + (m11*m33 - m13*m31)*(a12*(m22*m33 - m23*m32) - a22*(m12*m33 - m13*m32) + a32*(m12*m23 - m13*m22)) - (m21*m33 - m23*m31)*(a11*(m22*m33 - m23*m32) - a21*(m12*m33 - m13*m32) + a31*(m12*m23 - m13*m22)));
  result[2] = prefactor * ((m11*m22 - m12*m21)*(a13*(m22*m33 - m23*m32) - a23*(m12*m33 - m13*m32) + a33*(m12*m23 - m13*m22)) - (m11*m32 - m12*m31)*(a12*(m22*m33 - m23*m32) - a22*(m12*m33 - m13*m32) + a32*(m12*m23 - m13*m22)) + (m21*m32 - m22*m31)*(a11*(m22*m33 - m23*m32) - a21*(m12*m33 - m13*m32) + a31*(m12*m23 - m13*m22)));
  result[4] = prefactor * ((m11*m23 - m13*m21)*(a13*(m21*m33 - m23*m31) - a23*(m11*m33 - m13*m31) + a33*(m11*m23 - m13*m21)) - (m11*m33 - m13*m31)*(a12*(m21*m33 - m23*m31) - a22*(m11*m33 - m13*m31) + a32*(m11*m23 - m13*m21)) + (m21*m33 - m23*m31)*(a11*(m21*m33 - m23*m31) - a21*(m11*m33 - m13*m31) + a31*(m11*m23 - m13*m21)));
  result[5] = prefactor * (-(m11*m22 - m12*m21)*(a13*(m21*m33 - m23*m31) - a23*(m11*m33 - m13*m31) + a33*(m11*m23 - m13*m21)) + (m11*m32 - m12*m31)*(a12*(m21*m33 - m23*m31) - a22*(m11*m33 - m13*m31) + a32*(m11*m23 - m13*m21)) - (m21*m32 - m22*m31)*(a11*(m21*m33 - m23*m31) - a21*(m11*m33 - m13*m31) + a31*(m11*m23 - m13*m21)));
  result[8] = prefactor * ((m11*m22 - m12*m21)*(a13*(m21*m32 - m22*m31) - a23*(m11*m32 - m12*m31) + a33*(m11*m22 - m12*m21)) - (m11*m32 - m12*m31)*(a12*(m21*m32 - m22*m31) - a22*(m11*m32 - m12*m31) + a32*(m11*m22 - m12*m21)) + (m21*m32 - m22*m31)*(a11*(m21*m32 - m22*m31) - a21*(m11*m32 - m12*m31) + a31*(m11*m22 - m12*m21)));
  // we know that the diffusion tensor is symmetric. Therefore  J^-1 A J^-T  is symmetric, too.
  result[3] = result[1];
  result[6] = result[2];
  result[7] = result[5];
  // result[3] = prefactor * (-(m12*m23 - m13*m22)*(a13*(m21*m33 - m23*m31) - a23*(m11*m33 - m13*m31) + a33*(m11*m23 - m13*m21)) + (m12*m33 - m13*m32)*(a12*(m21*m33 - m23*m31) - a22*(m11*m33 - m13*m31) + a32*(m11*m23 - m13*m21)) - (m22*m33 - m23*m32)*(a11*(m21*m33 - m23*m31) - a21*(m11*m33 - m13*m31) + a31*(m11*m23 - m13*m21)));
  // result[6] = prefactor * ((m12*m23 - m13*m22)*(a13*(m21*m32 - m22*m31) - a23*(m11*m32 - m12*m31) + a33*(m11*m22 - m12*m21)) - (m12*m33 - m13*m32)*(a12*(m21*m32 - m22*m31) - a22*(m11*m32 - m12*m31) + a32*(m11*m22 - m12*m21)) + (m22*m33 - m23*m32)*(a11*(m21*m32 - m22*m31) - a21*(m11*m32 - m12*m31) + a31*(m11*m22 - m12*m21)));
  // result[7] = prefactor * (-(m11*m23 - m13*m21)*(a13*(m21*m32 - m22*m31) - a23*(m11*m32 - m12*m31) + a33*(m11*m22 - m12*m21)) + (m11*m33 - m13*m31)*(a12*(m21*m32 - m22*m31) - a22*(m11*m32 - m12*m31) + a32*(m11*m22 - m12*m21)) - (m21*m33 - m23*m31)*(a11*(m21*m32 - m22*m31) - a21*(m11*m32 - m12*m31) + a31*(m11*m22 - m12*m21)));

  return result;
}

// 1D integration factor
template<>
double computeIntegrationFactor<1>(const std::array<Vec3,1> &jacobian)
{
  return length<3>(jacobian[0]);
}

// 2D integration factor
template<>
double computeIntegrationFactor<2>(const std::array<Vec3,2> &jacobian)
{
  // rename input values
  const double m11 = jacobian[0][0];
  const double m21 = jacobian[0][1];
  const double m31 = jacobian[0][2];
  const double m12 = jacobian[1][0];
  const double m22 = jacobian[1][1];
  const double m32 = jacobian[1][2];
  return sqrt(sqr(m11*m22) + sqr(m11*m32) - 2*m11*m12*m21*m22 - 2*m11*m12*m31*m32
    + sqr(m12*m21) + sqr(m12*m31) + sqr(m21*m32) - 2*m21*m22*m31*m32 + sqr(m22*m31));
}

// 3D integration factor
template<>
double computeIntegrationFactor<3>(const std::array<Vec3,3> &jacobian)
{
  // rename input values
  const double m11 = jacobian[0][0];
  const double m21 = jacobian[0][1];
  const double m31 = jacobian[0][2];
  const double m12 = jacobian[1][0];
  const double m22 = jacobian[1][1];
  const double m32 = jacobian[1][2];
  const double m13 = jacobian[2][0];
  const double m23 = jacobian[2][1];
  const double m33 = jacobian[2][2];

  double determinant = m11*m22*m33 - m11*m23*m32 - m12*m21*m33 + m12*m23*m31 + m13*m21*m32 - m13*m22*m31;

  return fabs(determinant);
}

double applyTransformation(const std::array<double,9> &transformationMatrix, const Vec3 &vector1, const Vec3 &vector2)
{
  double result;

  // rename input values
  const double v11 = vector1[0];
  const double v12 = vector1[1];
  const double v13 = vector1[2];
  const double v21 = vector2[0];
  const double v22 = vector2[1];
  const double v23 = vector2[2];
  // 0 1 2
  // 3 4 5
  // 6 7 8
  const double m11 = transformationMatrix[0];
  const double m12 = transformationMatrix[1];
  const double m13 = transformationMatrix[2];
  const double m22 = transformationMatrix[4];
  const double m23 = transformationMatrix[5];
  const double m33 = transformationMatrix[8];

  //! computes v1^T * T * v2 where T is the symmetric transformation matrix
  //! computed by doc/compute_generalized_laplace.py

  // compute result
  result = v21*(m11*v11 + m12*v12 + m13*v13) + v22*(m12*v11 + m22*v12 + m23*v13) + v23*(m13*v11 + m23*v12 + m33*v13);
  return result;
}

double applyTransformation(const std::array<double,4> &transformationMatrix, const Vec2 &vector1, const Vec2 &vector2)
{
  double result;

  // rename input values
  const double v11 = vector1[0];
  const double v12 = vector1[1];
  const double v21 = vector2[0];
  const double v22 = vector2[1];
  // 0 1
  // 2 3
  const double m11 = transformationMatrix[0];
  const double m12 = transformationMatrix[1];
  const double m22 = transformationMatrix[3];

  //! computes v1^T * T * v2 where T is the symmetric transformation matrix
  //! computed by doc/compute_generalized_laplace.py

  // compute result
  result = v21*(m11*v11 + m12*v12) + v22*(m12*v11 + m22*v12);
  return result;
}

template<>
double computeDeterminant<3>(const Tensor2<3> &jacobian)
{
  // rename input values
  const double m11 = jacobian[0][0];
  const double m21 = jacobian[0][1];
  const double m31 = jacobian[0][2];
  const double m12 = jacobian[1][0];
  const double m22 = jacobian[1][1];
  const double m32 = jacobian[1][2];
  const double m13 = jacobian[2][0];
  const double m23 = jacobian[2][1];
  const double m33 = jacobian[2][2];

  return m11*m22*m33 - m11*m23*m32 - m12*m21*m33 + m12*m23*m31 + m13*m21*m32 - m13*m22*m31;
}

template<>
double computeDeterminant<2>(const Tensor2<2> &jacobian)
{
  // rename input values
  const double m11 = jacobian[0][0];
  const double m21 = jacobian[0][1];
  const double m12 = jacobian[1][0];
  const double m22 = jacobian[1][1];

  return m11*m22 - m12*m21;
}

template<>
Tensor2<3> computeSymmetricInverse<3>(const Tensor2<3> &matrix, double &determinant)
{
  // rename input values
  const double m11 = matrix[0][0];
  const double m21 = matrix[0][1];
  const double m31 = matrix[0][2];
  //const double m12 = matrix[1][0];
  const double m22 = matrix[1][1];
  const double m32 = matrix[1][2];
  const double m33 = matrix[2][2];

  determinant = m11*m22*m33 - m11*sqr(m32) - sqr(m21)*m33 + 2*m21*m31*m32 - m22*sqr(m31);
  double invDet = 1./determinant;

  std::array<Vec3,3> result;

  result[0][0] = invDet*(m22*m33 - sqr(m32));  // entry m11
  result[1][0] = invDet*(-m21*m33 + m31*m32);  // entry m12
  result[2][0] = invDet*(m21*m32 - m22*m31);   // entry m13
  result[0][1] = invDet*(-m21*m33 + m31*m32);  // entry m21
  result[1][1] = invDet*(m11*m33 - sqr(m31));  // entry m22
  result[2][1] = invDet*(-m11*m32 + m21*m31);  // entry m23
  result[0][2] = invDet*(m21*m32 - m22*m31);   // entry m31
  result[1][2] = invDet*(-m11*m32 + m21*m31);  // entry m32
  result[2][2] = invDet*(m11*m22 - sqr(m21));  // entry m33

  return result;
}

template<>
Tensor2<2> computeSymmetricInverse<2>(const Tensor2<2> &matrix, double &determinant)
{
  // rename input values
  const double m11 = matrix[0][0];
  const double m21 = matrix[0][1];
  //const double m12 = matrix[1][0];
  const double m22 = matrix[1][1];

  determinant = m11*m22 - sqr(m21);
  double invDet = 1./determinant;

  Tensor2<2> result;

  result[0][0] = invDet * m22;  // entry m11
  result[1][0] = invDet * -m21;  // entry m12
  result[0][1] = invDet * -m21;  // entry m21
  result[1][1] = invDet * m11;  // entry m22

  return result;
}

template<>
Tensor2<3> computeInverse<3>(const Tensor2<3> &matrix, double &determinant)
{
  // matrices are stored column-major

  // rename input values
  const double m11 = matrix[0][0];
  const double m21 = matrix[0][1];
  const double m31 = matrix[0][2];
  const double m12 = matrix[1][0];
  const double m22 = matrix[1][1];
  const double m32 = matrix[1][2];
  const double m13 = matrix[2][0];
  const double m23 = matrix[2][1];
  const double m33 = matrix[2][2];

  determinant =  m11*m22*m33 - m11*m23*m32 - m12*m21*m33 + m12*m23*m31 + m13*m21*m32 - m13*m22*m31;
  double invDet = 1./determinant;

  std::array<Vec3,3> result;

  result[0][0] = invDet*(m22*m33 - m23*m32);   // entry m11
  result[1][0] = invDet*(-m12*m33 + m13*m32);  // entry m12
  result[2][0] = invDet*(m12*m23 - m13*m22);   // entry m13
  result[0][1] = invDet*(-m21*m33 + m23*m31);  // entry m21
  result[1][1] = invDet*(m11*m33 - m13*m31);   // entry m22
  result[2][1] = invDet*(-m11*m23 + m13*m21);  // entry m23
  result[0][2] = invDet*(m21*m32 - m22*m31);   // entry m31
  result[1][2] = invDet*(-m11*m32 + m12*m31);  // entry m32
  result[2][2] = invDet*(m11*m22 - m12*m21);   // entry m33

  return result;
}

template<>
Tensor2<2> computeInverse<2>(const Tensor2<2> &matrix, double &determinant)
{
  // matrices are stored column-major

  // rename input values
  const double m11 = matrix[0][0];
  const double m21 = matrix[0][1];
  const double m12 = matrix[1][0];
  const double m22 = matrix[1][1];

  determinant =  m11*m22 - m12*m21;
  double invDet = 1./determinant;

  Tensor2<2> result;

  result[0][0] = invDet*(m22);   // entry m11
  result[1][0] = invDet*(-m12);  // entry m12
  result[0][1] = invDet*(-m21);  // entry m21
  result[1][1] = invDet*(m11);   // entry m22

  return result;
}

template<>
Tensor2<1> computeInverse<1>(const Tensor2<1> &matrix, double &determinant)
{
  determinant = matrix[0][0];
  
  Tensor2<1> result;
  result[0][0] = 1./determinant;
  return result;
}

template<>
Tensor2<3> computeCofactorMatrix<3>(const Tensor2<3> &matrix)
{
  // matrices are stored column-major

  // rename input values
  const double m11 = matrix[0][0];
  const double m21 = matrix[0][1];
  const double m31 = matrix[0][2];
  const double m12 = matrix[1][0];
  const double m22 = matrix[1][1];
  const double m32 = matrix[1][2];
  const double m13 = matrix[2][0];
  const double m23 = matrix[2][1];
  const double m33 = matrix[2][2];

  std::array<Vec3,3> result;

  result[0][0] = m22*m33 - m23*m32;   // entry m11
  result[1][0] = -m21*m33 + m23*m31;  // entry m12
  result[2][0] = m21*m32 - m22*m31;   // entry m13
  result[0][1] = -m12*m33 + m13*m32;  // entry m21
  result[1][1] = m11*m33 - m13*m31;   // entry m22
  result[2][1] = -m11*m32 + m12*m31;  // entry m23
  result[0][2] = m12*m23 - m13*m22;   // entry m31
  result[1][2] = -m11*m23 + m13*m21;  // entry m32
  result[2][2] = m11*m22 - m12*m21;   // entry m33

  return result;
}

template<>
Tensor2<2> computeCofactorMatrix<2>(const Tensor2<2> &matrix)
{
  // matrices are stored column-major

  // rename input values
  const double m11 = matrix[0][0];
  const double m21 = matrix[0][1];
  const double m12 = matrix[1][0];
  const double m22 = matrix[1][1];

  Tensor2<2> result;

  result[0][0] = m22;   // entry m11
  result[1][0] = -m21;  // entry m12
  result[0][1] = -m12;  // entry m21
  result[1][1] = m11;   // entry m22

  return result;
}

//! transform a 3xD2 matrix to a DxD matrix by filling up with identity entries
template<>
std::array<std::array<double,3>,3> transformToDxD<3,3>(const std::array<Vec3,3> &matrix)
{
  return matrix;
}

template<>
std::array<std::array<double,3>,3> transformToDxD<3,2>(const std::array<Vec3,2> &matrix)
{
  return std::array<Vec3,3>({matrix[0], matrix[1], Vec3({0.0, 0.0, 1.0})});
}

template<>
std::array<std::array<double,2>,2> transformToDxD<2,2>(const std::array<Vec3,2> &matrix)
{
  return std::array<std::array<double,2>,2>({
    std::array<double,2>({matrix[0][0], matrix[0][1]}),
    std::array<double,2>({matrix[1][0], matrix[1][1]})});
}

template<>
std::array<std::array<double,3>,3> transformToDxD<3,1>(const std::array<Vec3,1> &matrix)
{
  return std::array<Vec3,3>({matrix[0], Vec3({0.0, 1.0, 0.0}), Vec3({0.0, 0.0, 1.0})});
}

template<>
std::array<std::array<double,1>,1> transformToDxD<1,1>(const std::array<Vec3,1> &matrix)
{
  return std::array<std::array<double,1>,1>({std::array<double,1>({matrix[0][0]})});
}

template<>
VecD<3> transformToD<3,3>(const VecD<3> &vector)
{
  return vector;
}

template<>
VecD<3> transformToD<3,2>(const VecD<2> &vector)
{
  return VecD<3>({vector[0], vector[1], 0.0});
}

template<>
VecD<2> transformToD<2,3>(const VecD<3> &vector)
{
  return VecD<2>({vector[0], vector[1]});
}

template<>
VecD<1> transformToD<1,3>(const VecD<3> &vector)
{
  return VecD<1>({vector[0]});
}

//! compute 3D cross product
Vec3 cross(const Vec3 &vector1, const Vec3 &vector2)
{
  return Vec3{
    vector1[1]*vector2[2] - vector1[2]*vector2[1], 
    vector1[2]*vector2[0] - vector1[0]*vector2[2], 
    vector1[0]*vector2[1] - vector1[1]*vector2[0]};
}

//! compute dot product
double dot(const Vec3 &vector1, const Vec3 &vector2)
{
  return vector1[0]*vector2[0] + vector1[1]*vector2[1] + vector1[2]*vector2[2];
}

bool isSubsequenceOf(std::vector<int> a, std::vector<int> b, size_t &subsequenceAStartPos)
{
  if (b.empty())
    return true;

  // find the matching entry in vector a
  bool matchFound = false;
  size_t aIndex=0;
  for (; aIndex<a.size(); aIndex++)
  {
    if (a[aIndex] == b[0])
    {
      subsequenceAStartPos = aIndex;
      matchFound = true;
      break;
    }
  }

  if (!matchFound)
    return false;

  for (size_t bIndex=1; bIndex<b.size(); bIndex++)
  {
    aIndex++;
    if (aIndex >= a.size())
      return false;

    if (a[aIndex] != b[bIndex])
      return false;
  }
  return true;
}

int permutation(int i, int j, int k)
{
  if ((i==1 && j==2 && k==3) || (i==2 && j==3 && k==1) || (i==3 && j==1 && k==2))
    return 1;
  else if ((i==3 && j==2 && k==1) || (i==1 && j==3 && k==2) || (i==2 && j==1 && k==3))
    return -1;
  return 0;
}

void quadrilateralGetPointCoordinates(const std::array<Vec3,4> geometryValues, const Vec3 point, Vec2 &xi)
{
  // derivation using sympy in script invert_mapping.py
  const double xp1 = point[0];
  const double xp2 = point[1];

  double x11 = geometryValues[0][0];
  double x12 = geometryValues[0][1];

  const double x21 = geometryValues[1][0];
  const double x22 = geometryValues[1][1];

  const double x31 = geometryValues[2][0];
  const double x32 = geometryValues[2][1];

  const double x41 = geometryValues[3][0];
  const double x42 = geometryValues[3][1];

  // compute analytic solution for xi
  // avoid division by 0
  const double divisor = (x11*x32 - x11*x42 - x12*x31 + x12*x41 - x21*x32 + x21*x42 + x22*x31 - x22*x41);
  const double eps = 1e-12;
  if (fabs(divisor) < eps)
  {
    x11 += eps;
    x12 += 0.8*eps;
  }

  double xi1 = 0.5*(2*x11*x32 - x11*x42 - 2*x12*x31 + x12*x41 - x21*x32 + x22*x31 + xp1*(x12 - x22 - x32 + x42) + xp2*(-x11 + x21 + x31 - x41) + std::sqrt(std::pow(x11, 2)*std::pow(x42, 2) - 2*std::pow(x11, 2)*x42*xp2 + std::pow(x11, 2)*std::pow(xp2, 2) - 2*x11*x12*x41*x42 + 2*x11*x12*x41*xp2 + 2*x11*x12*x42*xp1 - 2*x11*x12*xp1*xp2 - 2*x11*x21*x32*x42 + 2*x11*x21*x32*xp2 + 2*x11*x21*x42*xp2 - 2*x11*x21*std::pow(xp2, 2) - 2*x11*x22*x31*x42 + 2*x11*x22*x31*xp2 + 4*x11*x22*x32*x41 - 4*x11*x22*x32*xp1 - 4*x11*x22*x41*xp2 + 2*x11*x22*x42*xp1 + 2*x11*x22*xp1*xp2 + 2*x11*x31*x42*xp2 - 2*x11*x31*std::pow(xp2, 2) - 4*x11*x32*x41*xp2 + 2*x11*x32*x42*xp1 + 2*x11*x32*xp1*xp2 + 2*x11*x41*x42*xp2 + 2*x11*x41*std::pow(xp2, 2) - 2*x11*std::pow(x42, 2)*xp1 - 2*x11*x42*xp1*xp2 + std::pow(x12, 2)*std::pow(x41, 2) - 2*std::pow(x12, 2)*x41*xp1 + std::pow(x12, 2)*std::pow(xp1, 2) + 4*x12*x21*x31*x42 - 4*x12*x21*x31*xp2 - 2*x12*x21*x32*x41 + 2*x12*x21*x32*xp1 + 2*x12*x21*x41*xp2 - 4*x12*x21*x42*xp1 + 2*x12*x21*xp1*xp2 - 2*x12*x22*x31*x41 + 2*x12*x22*x31*xp1 + 2*x12*x22*x41*xp1 - 2*x12*x22*std::pow(xp1, 2) + 2*x12*x31*x41*xp2 - 4*x12*x31*x42*xp1 + 2*x12*x31*xp1*xp2 + 2*x12*x32*x41*xp1 - 2*x12*x32*std::pow(xp1, 2) - 2*x12*std::pow(x41, 2)*xp2 + 2*x12*x41*x42*xp1 - 2*x12*x41*xp1*xp2 + 2*x12*x42*std::pow(xp1, 2) + std::pow(x21, 2)*std::pow(x32, 2) - 2*std::pow(x21, 2)*x32*xp2 + std::pow(x21, 2)*std::pow(xp2, 2) - 2*x21*x22*x31*x32 + 2*x21*x22*x31*xp2 + 2*x21*x22*x32*xp1 - 2*x21*x22*xp1*xp2 + 2*x21*x31*x32*xp2 - 4*x21*x31*x42*xp2 + 2*x21*x31*std::pow(xp2, 2) - 2*x21*std::pow(x32, 2)*xp1 + 2*x21*x32*x41*xp2 + 2*x21*x32*x42*xp1 - 2*x21*x32*xp1*xp2 - 2*x21*x41*std::pow(xp2, 2) + 2*x21*x42*xp1*xp2 + std::pow(x22, 2)*std::pow(x31, 2) - 2*std::pow(x22, 2)*x31*xp1 + std::pow(x22, 2)*std::pow(xp1, 2) - 2*x22*std::pow(x31, 2)*xp2 + 2*x22*x31*x32*xp1 + 2*x22*x31*x41*xp2 + 2*x22*x31*x42*xp1 - 2*x22*x31*xp1*xp2 - 4*x22*x32*x41*xp1 + 2*x22*x32*std::pow(xp1, 2) + 2*x22*x41*xp1*xp2 - 2*x22*x42*std::pow(xp1, 2) + std::pow(x31, 2)*std::pow(xp2, 2) - 2*x31*x32*xp1*xp2 - 2*x31*x41*std::pow(xp2, 2) + 2*x31*x42*xp1*xp2 + std::pow(x32, 2)*std::pow(xp1, 2) + 2*x32*x41*xp1*xp2 - 2*x32*x42*std::pow(xp1, 2) + std::pow(x41, 2)*std::pow(xp2, 2) - 2*x41*x42*xp1*xp2 + std::pow(x42, 2)*std::pow(xp1, 2)))/(x11*x32 - x11*x42 - x12*x31 + x12*x41 - x21*x32 + x21*x42 + x22*x31 - x22*x41);


  // avoid division by 0
  const double divisor2 = (x11*xi1 - x11 - x21*xi1 - x31*xi1 + x31 + x41*xi1);
  if (fabs(divisor2) < eps)
  {
    x11 -= 0.6*eps;
    xi1 += 0.8*eps;
  }
  const double xi2 = (x11*xi1 - x11 - x21*xi1 + xp1)/(x11*xi1 - x11 - x21*xi1 - x31*xi1 + x31 + x41*xi1);

  xi[0] = xi1;
  xi[1] = xi2;
}

template<>
void rotateMatrix<1>(Matrix<1,1> &matrix, Vec3 directionVector)
{
  // 1D case does not make sense
  assert(false);
}

template<>
void rotateMatrix<2>(Matrix<2,2> &matrix, Vec3 directionVector)
{
  // derivation in compute_rotation_tensor.py
  const double b1 = directionVector[0];
  const double b2 = directionVector[1];

  Matrix<2,2> rotationMatrix(
  {
    b1, -b2,
    b2, -b1
  });

  const double determinant = (b1 - b2)*(b1 + b2);  // b1^2 - b2^2
  Matrix<2,2> rotationMatrixInverse(
  {
    b1/determinant, -b2/determinant,
    b2/determinant, -b1/determinant
  });

  matrix = rotationMatrixInverse * matrix * rotationMatrix;
}

template<>
void rotateMatrix<3>(Matrix<3,3> &matrix, Vec3 directionVector)
{
  const double b1 = directionVector[0];
  const double b2 = directionVector[1];
  const double b3 = directionVector[2];

  Matrix<3,3> basis(
  {
    // normalized(b)                       normalized(e_3 x b)                        normalized(b x e_3 x b)
    b1/sqrt(sqr(b1) + sqr(b2) + sqr(b3)),  -b2/sqrt(sqr(b1) + sqr(b2)),               -b1*b3/sqrt(sqr(b1)*sqr(b3) + sqr(b2)*sqr(b3) + sqr(sqr(b1) + sqr(b2))),
    b2/sqrt(sqr(b1) + sqr(b2) + sqr(b3)),   b1/sqrt(sqr(b1) + sqr(b2)),               -b2*b3/sqrt(sqr(b1)*sqr(b3) + sqr(b2)*sqr(b3) + sqr(sqr(b1) + sqr(b2))),
    b3/sqrt(sqr(b1) + sqr(b2) + sqr(b3)),                            0,  (sqr(b1) + sqr(b2))/sqrt(sqr(b1)*sqr(b3) + sqr(b2)*sqr(b3) + sqr(sqr(b1) + sqr(b2)))
  });

  Matrix<3,3> basis_T(
  {
    basis(0,0), basis(1,0), basis(2,0),
    basis(0,1), basis(1,1), basis(2,1),
    basis(0,2), basis(1,2), basis(2,2)
  });

  /*VLOG(1) << "direction: " << directionVector;
  VLOG(1) << "rotation basis: " << basis;

  VLOG(1) << "matrix before rotation: " << matrix;
  VLOG(1) << "   basis * matrix * basis_T: " << basis * matrix * basis_T;
  */

  matrix = basis * matrix * basis_T;

  //VLOG(1) << "matrix after rotation: " << matrix;
}

double estimateMaximumEigenvalue(const Tensor2<3> &matrix)
{
  Vec3 v({1.0,0.0,0.0});

  double normVPrevious = 0;
  double normV = 1;

  // power iteration
  for (int i = 0; i < 15 && fabs(normVPrevious - normV) >= 1e-5; i++)
  {
    normVPrevious = normV;

    v = matrix*v;

    // normalize vector
    normV = norm<3>(v);
    v /= normV;

    VLOG(1) << " eigenvalue estimation, it " << i << ", v: " << v << ", value: " << normV << ", error: " << fabs(normVPrevious - normV);
  }

  return normV;
}

double estimateConditionNumber(const Tensor2<3> &matrix, const Tensor2<3> &inverseMatrix)
{
  double maximumEigenvalue = estimateMaximumEigenvalue(matrix);
  double minimumEigenvalue = 1./estimateMaximumEigenvalue(inverseMatrix);

  return maximumEigenvalue / minimumEigenvalue;
}

}  // namespace
