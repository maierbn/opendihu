#include "utility/math_utility.h"

#include <cmath>

namespace MathUtility
{

const double INVERSE_REGULARIZATION_TOLERANCE = 1e-10;
const double INVERSE_REGULARIZATION_EPSILON = 1e-10;


template<typename double_v_t>
double_v_t length(const VecD<1,double_v_t> node)
{
  return node[0];
}

template<typename double_v_t>
double_v_t length(const VecD<2,double_v_t> node)
{
  return sqrt(sqr(node[0]) + sqr(node[1]));
}

template<typename double_v_t>
double_v_t length(const VecD<3,double_v_t> node)
{
  return sqrt(sqr(node[0]) + sqr(node[1]) + sqr(node[2]));
}

template<int D, typename double_v_t>
double_v_t norm(const VecD<D,double_v_t> node)
{
  return length(node);
}

template<int D>
VecD<D> normalized(VecD<D> &vector)
{
  double factor = 1./norm<D>(vector);
  return vector * factor;
}

template<int D>
void normalize(VecD<D> &vector)
{
  vector *= 1./norm<D>(vector);
}

template<typename T>
void readPoint(T &file, Vec3 &point)
{
  union
  {
    char c[8];
    double d;
  };
  for (int i = 0; i < 3; i++)
  {
    file.read(c, 8);
    point[i] = d;
  }
}

template<typename T>
void writePoint(T &file, Vec3 &point)
{
  union
  {
    char c[8];
    double d;
  };
  for (int i = 0; i < 3; i++)
  {
    d = point[i];
    file.write(c, 8);
  }
}

template<int D>
double distance(const VecD<D> node1, const VecD<D> node2)
{
  double result = 0;
  for (int componentNo = 0; componentNo < D; componentNo++)
  {
    result += sqr(node1[componentNo]-node2[componentNo]);
  }
  return sqrt(result);
}

// 1D integration factor
template<typename double_v_t>
double_v_t computeIntegrationFactor(const std::array<VecD<3,double_v_t>,1> &jacobian)
{
  LOG(DEBUG) << "1D integration factor, jacobian[0]: " << jacobian;
  return length(jacobian[0]);
}

// 2D integration factor
template<typename double_v_t>
double_v_t computeIntegrationFactor(const std::array<VecD<3,double_v_t>,2> &jacobian)
{
  // rename input values
  const double_v_t m11 = jacobian[0][0];
  const double_v_t m21 = jacobian[0][1];
  const double_v_t m31 = jacobian[0][2];
  const double_v_t m12 = jacobian[1][0];
  const double_v_t m22 = jacobian[1][1];
  const double_v_t m32 = jacobian[1][2];

  /*
  LOG(DEBUG) << "2D integration factor, m11:" << m11 << ", m21: " << m21 << ", m32: " << m32 << ", sqr(m11*m22): " << sqr(m11*m22) << ", sqrt "  << sqr(m11*m22) + sqr(m11*m32) - 2*m11*m12*m21*m22 - 2*m11*m12*m31*m32
    + sqr(m12*m21) + sqr(m12*m31) + sqr(m21*m32) - 2*m21*m22*m31*m32 + sqr(m22*m31) << ", result: " << sqrt(sqr(m11*m22) + sqr(m11*m32) - 2*m11*m12*m21*m22 - 2*m11*m12*m31*m32
    + sqr(m12*m21) + sqr(m12*m31) + sqr(m21*m32) - 2*m21*m22*m31*m32 + sqr(m22*m31));
    */

  return sqrt(sqr(m11*m22) + sqr(m11*m32) - 2*m11*m12*m21*m22 - 2*m11*m12*m31*m32
    + sqr(m12*m21) + sqr(m12*m31) + sqr(m21*m32) - 2*m21*m22*m31*m32 + sqr(m22*m31));
}

// 3D integration factor
template<typename double_v_t>
double_v_t computeIntegrationFactor(const std::array<VecD<3,double_v_t>,3> &jacobian)
{
  // rename input values
  const double_v_t m11 = jacobian[0][0];
  const double_v_t m21 = jacobian[0][1];
  const double_v_t m31 = jacobian[0][2];
  const double_v_t m12 = jacobian[1][0];
  const double_v_t m22 = jacobian[1][1];
  const double_v_t m32 = jacobian[1][2];
  const double_v_t m13 = jacobian[2][0];
  const double_v_t m23 = jacobian[2][1];
  const double_v_t m33 = jacobian[2][2];

  double_v_t determinant = m11*m22*m33 - m11*m23*m32 - m12*m21*m33 + m12*m23*m31 + m13*m21*m32 - m13*m22*m31;

  return abs(determinant);
}

// 3D matrix determinant
template<typename double_v_t>
double_v_t computeDeterminant(const Tensor2<3,double_v_t> &jacobian)
{
  // rename input values
  const double_v_t m11 = jacobian[0][0];
  const double_v_t m21 = jacobian[0][1];
  const double_v_t m31 = jacobian[0][2];
  const double_v_t m12 = jacobian[1][0];
  const double_v_t m22 = jacobian[1][1];
  const double_v_t m32 = jacobian[1][2];
  const double_v_t m13 = jacobian[2][0];
  const double_v_t m23 = jacobian[2][1];
  const double_v_t m33 = jacobian[2][2];

  return m11*m22*m33 - m11*m23*m32 - m12*m21*m33 + m12*m23*m31 + m13*m21*m32 - m13*m22*m31;
}

// 2D matrix determinant
template<typename double_v_t>
double_v_t computeDeterminant(const Tensor2<2> &jacobian)
{
  // rename input values
  const double_v_t m11 = jacobian[0][0];
  const double_v_t m21 = jacobian[0][1];
  const double_v_t m12 = jacobian[1][0];
  const double_v_t m22 = jacobian[1][1];

  return m11*m22 - m12*m21;
}

// 3D symmetric matrix inverse
template<typename double_v_t>
Tensor2<3,double_v_t> computeSymmetricInverse(const Tensor2<3,double_v_t> &matrix, double_v_t &determinant)
{
  // rename input values
        double_v_t m11 = matrix[0][0];
  const double_v_t m21 = matrix[0][1];
  const double_v_t m31 = matrix[0][2];
  //const double_v_t m12 = matrix[1][0];
        double_v_t m22 = matrix[1][1];
  const double_v_t m32 = matrix[1][2];
        double_v_t m33 = matrix[2][2];

  determinant = m11*m22*m33 - m11*sqr(m32) - sqr(m21)*m33 + 2*m21*m31*m32 - m22*sqr(m31);

  // regularize matrix if near singular, by adding ε*I (small values on diagonal)
  double_v_t regularizationEpsilon = 0;
  // where operator is documented here: https://vcdevel.github.io/Vc-1.4.1/group__Utilities.html#gaa18ac68167ac7614731134de7364a1d5
  Vc::where(abs(determinant) < INVERSE_REGULARIZATION_TOLERANCE) | regularizationEpsilon = INVERSE_REGULARIZATION_EPSILON;
  m11 += regularizationEpsilon;
  m22 += regularizationEpsilon;
  m33 += regularizationEpsilon;

  determinant = m11*m22*m33 - m11*sqr(m32) - sqr(m21)*m33 + 2*m21*m31*m32 - m22*sqr(m31);
  double_v_t invDet = 1./determinant;

  Tensor2<3,double_v_t> result;

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

// 2D symmetric matrix invers
template<typename double_v_t>
Tensor2<2,double_v_t> computeSymmetricInverse(const Tensor2<2,double_v_t> &matrix, double_v_t &determinant)
{
  // rename input values
        double_v_t m11 = matrix[0][0];
  const double_v_t m21 = matrix[0][1];
  //const double_v_t m12 = matrix[1][0];
        double_v_t m22 = matrix[1][1];

  determinant = m11*m22 - sqr(m21);

  // regularize matrix if near singular, by adding ε*I (small values on diagonal)
  double_v_t regularizationEpsilon = 0;
  // where operator is documented here: https://vcdevel.github.io/Vc-1.4.1/group__Utilities.html#gaa18ac68167ac7614731134de7364a1d5
  Vc::where(abs(determinant) < INVERSE_REGULARIZATION_TOLERANCE) | regularizationEpsilon = INVERSE_REGULARIZATION_EPSILON;

  m11 += regularizationEpsilon;
  m22 += regularizationEpsilon;

  determinant = m11*m22 - sqr(m21);
  double_v_t invDet = 1./determinant;

  Tensor2<2,double_v_t> result;

  result[0][0] = invDet * m22;  // entry m11
  result[1][0] = invDet * -m21;  // entry m12
  result[0][1] = invDet * -m21;  // entry m21
  result[1][1] = invDet * m11;  // entry m22

  return result;
}

// 1D inverse
template<typename double_v_t>
Tensor2<1,double_v_t> computeInverse(const Tensor2<1,double_v_t> &matrix, double_v_t &determinant)
{
  // regularize matrix if near singular, by adding ε*I (small values on diagonal)
  double_v_t regularizationEpsilon = 0;
  // where operator is documented here: https://vcdevel.github.io/Vc-1.4.1/group__Utilities.html#gaa18ac68167ac7614731134de7364a1d5
  Vc::where(abs(matrix[0][0]) < INVERSE_REGULARIZATION_TOLERANCE) | regularizationEpsilon = INVERSE_REGULARIZATION_EPSILON;

  determinant = matrix[0][0] + regularizationEpsilon;

  Tensor2<1,double_v_t> result;
  result[0][0] = 1./determinant;
  return result;
}

// 2D inverse
template<typename double_v_t>
Tensor2<2,double_v_t> computeInverse(const Tensor2<2,double_v_t> &matrix, double_v_t &determinant)
{
  // matrices are stored column-major

  // rename input values
        double_v_t m11 = matrix[0][0];
  const double_v_t m21 = matrix[0][1];
  const double_v_t m12 = matrix[1][0];
        double_v_t m22 = matrix[1][1];

  determinant =  m11*m22 - m12*m21;

  // regularize matrix if near singular, by adding ε*I (small values on diagonal)
  double_v_t regularizationEpsilon = 0;
  // where operator is documented here: https://vcdevel.github.io/Vc-1.4.1/group__Utilities.html#gaa18ac68167ac7614731134de7364a1d5
  Vc::where(abs(determinant) < INVERSE_REGULARIZATION_TOLERANCE) | regularizationEpsilon = INVERSE_REGULARIZATION_EPSILON;

  m11 += regularizationEpsilon;
  m22 += regularizationEpsilon;

  determinant =  m11*m22 - m12*m21;
  double_v_t invDet = 1./determinant;

  Tensor2<2,double_v_t> result;

  result[0][0] = invDet*(m22);   // entry m11
  result[1][0] = invDet*(-m12);  // entry m12
  result[0][1] = invDet*(-m21);  // entry m21
  result[1][1] = invDet*(m11);   // entry m22

  return result;
}

// 3D inverse
template<typename double_v_t>
Tensor2<3,double_v_t> computeInverse(const Tensor2<3,double_v_t> &matrix, double_v_t &determinant)
{
  // rename input values
  // matrices are stored column-major
        double_v_t m11 = matrix[0][0];
  const double_v_t m21 = matrix[0][1];
  const double_v_t m31 = matrix[0][2];
  const double_v_t m12 = matrix[1][0];
        double_v_t m22 = matrix[1][1];
  const double_v_t m32 = matrix[1][2];
  const double_v_t m13 = matrix[2][0];
  const double_v_t m23 = matrix[2][1];
        double_v_t m33 = matrix[2][2];

  determinant =  m11*m22*m33 - m11*m23*m32 - m12*m21*m33 + m12*m23*m31 + m13*m21*m32 - m13*m22*m31;

  // regularize matrix if near singular, by adding ε*I (small values on diagonal)
  double_v_t regularizationEpsilon = 0;

  // where operator is documented here: https://vcdevel.github.io/Vc-1.4.1/group__Utilities.html#gaa18ac68167ac7614731134de7364a1d5
  Vc::where(abs(determinant) < INVERSE_REGULARIZATION_TOLERANCE) | regularizationEpsilon = INVERSE_REGULARIZATION_EPSILON;

  m11 += regularizationEpsilon;
  m22 += regularizationEpsilon;
  m33 += regularizationEpsilon;

  determinant =  m11*m22*m33 - m11*m23*m32 - m12*m21*m33 + m12*m23*m31 + m13*m21*m32 - m13*m22*m31;
  double_v_t invDet = 1./determinant;

  Tensor2<3,double_v_t> result;

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

template<typename double_v_t=double>
std::array<double_v_t,9> computeTransformationMatrixAndDeterminant(const std::array<VecD<3,double_v_t>,3> &jacobian, double_v_t &determinant)
{
  // rename input values
  const double_v_t m11 = jacobian[0][0];
  const double_v_t m21 = jacobian[0][1];
  const double_v_t m31 = jacobian[0][2];
  const double_v_t m12 = jacobian[1][0];
  const double_v_t m22 = jacobian[1][1];
  const double_v_t m32 = jacobian[1][2];
  const double_v_t m13 = jacobian[2][0];
  const double_v_t m23 = jacobian[2][1];
  const double_v_t m33 = jacobian[2][2];

  // the following code is copy&pasted from the output of ./invert_mapping.py
  determinant = m11*m22*m33 - m11*m23*m32 - m12*m21*m33 + m12*m23*m31 + m13*m21*m32 - m13*m22*m31;
  double_v_t prefactor = 1./sqr(determinant);

  // 0 1 2
  // 3 4 5
  // 6 7 8

  // compute matrix entries of symmetric matrix
  std::array<double_v_t, 9> result;
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

template<typename double_v_t=double>
std::array<double_v_t,9> computeTransformationDiffusionMatrixAndDeterminant(const std::array<VecD<3,double_v_t>,3> &jacobian, const Matrix<3,3,double_v_t> &diffusionTensor, double_v_t &determinant)
{
  // rename input values
  const double_v_t m11 = jacobian[0][0];
  const double_v_t m21 = jacobian[0][1];
  const double_v_t m31 = jacobian[0][2];
  const double_v_t m12 = jacobian[1][0];
  const double_v_t m22 = jacobian[1][1];
  const double_v_t m32 = jacobian[1][2];
  const double_v_t m13 = jacobian[2][0];
  const double_v_t m23 = jacobian[2][1];
  const double_v_t m33 = jacobian[2][2];

  const double_v_t a11 = diffusionTensor(0,0);
  const double_v_t a21 = diffusionTensor(0,1);
  const double_v_t a31 = diffusionTensor(0,2);
  const double_v_t a12 = diffusionTensor(1,0);
  const double_v_t a22 = diffusionTensor(1,1);
  const double_v_t a32 = diffusionTensor(1,2);
  const double_v_t a13 = diffusionTensor(2,0);
  const double_v_t a23 = diffusionTensor(2,1);
  const double_v_t a33 = diffusionTensor(2,2);

  // the following code is copy&pasted from the output of ./invert_mapping.py
  determinant = m11*m22*m33 - m11*m23*m32 - m12*m21*m33 + m12*m23*m31 + m13*m21*m32 - m13*m22*m31;
  double_v_t prefactor = 1./sqr(determinant);

  // 0 1 2
  // 3 4 5
  // 6 7 8

  // compute matrix entries of symmetric matrix  J^-1 A J^-T  assuming that A is symmetric
  std::array<double_v_t, 9> result;
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

template<typename double_v1_t,typename double_v2_t,typename double_v3_t>
double_v1_t applyTransformation(const std::array<double_v1_t,9> &transformationMatrix, const VecD<3,double_v2_t> &vector1, const VecD<3,double_v3_t> &vector2)
{
  double_v1_t result;

  // rename input values
  const double_v2_t v11 = vector1[0];
  const double_v2_t v12 = vector1[1];
  const double_v2_t v13 = vector1[2];
  const double_v3_t v21 = vector2[0];
  const double_v3_t v22 = vector2[1];
  const double_v3_t v23 = vector2[2];
  // 0 1 2
  // 3 4 5
  // 6 7 8
  const double_v1_t m11 = transformationMatrix[0];
  const double_v1_t m12 = transformationMatrix[1];
  const double_v1_t m13 = transformationMatrix[2];
  const double_v1_t m22 = transformationMatrix[4];
  const double_v1_t m23 = transformationMatrix[5];
  const double_v1_t m33 = transformationMatrix[8];

  //! computes v1^T * T * v2 where T is the symmetric transformation matrix
  //! computed by doc/compute_generalized_laplace.py

  // compute result
  result = v21*(m11*v11 + m12*v12 + m13*v13) + v22*(m12*v11 + m22*v12 + m23*v13) + v23*(m13*v11 + m23*v12 + m33*v13);
  return result;
}

template<typename double_v1_t,typename double_v2_t,typename double_v3_t>
double_v1_t applyTransformation(const std::array<double_v1_t,4> &transformationMatrix, const VecD<2,double_v2_t> &vector1, const VecD<2,double_v3_t> &vector2)
{
  double_v1_t result;

  // rename input values
  const double_v2_t v11 = vector1[0];
  const double_v2_t v12 = vector1[1];
  const double_v3_t v21 = vector2[0];
  const double_v3_t v22 = vector2[1];
  // 0 1
  // 2 3
  const double_v1_t m11 = transformationMatrix[0];
  const double_v1_t m12 = transformationMatrix[1];
  const double_v1_t m22 = transformationMatrix[3];

  //! computes v1^T * T * v2 where T is the symmetric transformation matrix
  //! computed by doc/compute_generalized_laplace.py

  // compute result
  result = v21*(m11*v11 + m12*v12) + v22*(m12*v11 + m22*v12);
  return result;
}

template<typename double_v_t>
void rotateMatrix(Matrix<2,2,double_v_t> &matrix, VecD<2,double_v_t> directionVector)
{
  // derivation in compute_rotation_tensor.py
  const double_v_t b1 = directionVector[0];
  const double_v_t b2 = directionVector[1];

  Matrix<2,2,double_v_t> rotationMatrix(
  {
    b1, -b2,
    b2, -b1
  });

  const double_v_t determinant = (b1 - b2)*(b1 + b2);  // b1^2 - b2^2
  Matrix<2,2,double_v_t> rotationMatrixInverse(
  {
    b1/determinant, -b2/determinant,
    b2/determinant, -b1/determinant
  });

  matrix = rotationMatrixInverse * matrix * rotationMatrix;
}

template<typename double_v_t>
void rotateMatrix(Matrix<3,3,double_v_t> &matrix, VecD<3,double_v_t> directionVector)
{
  const double_v_t b1 = directionVector[0];
  const double_v_t b2 = directionVector[1];
  const double_v_t b3 = directionVector[2];

  const double_v_t bNorm = sqrt(sqr(b1) + sqr(b2) + sqr(b3));     // |b|
  const double_v_t e3bNorm = sqrt(sqr(b1) + sqr(b2));             // |e_3 x b|
  const double_v_t be3bNorm = sqrt(sqr(b1)*sqr(b3) + sqr(b2)*sqr(b3) + sqr(sqr(b1) + sqr(b2)));   // |b x (e_3 x b)|

  if (MathUtility::abs(bNorm) < 1e-15)
  {
    LOG(ERROR) << "Trying to transform a matrix to direction vector " << directionVector;
    return;
  }

  // if one of the new basis vectors is zero (this happens, if directionVector == [0,0,1]), use e2 instead of e3
  if (MathUtility::abs(e3bNorm) < 1e-15 || MathUtility::abs(be3bNorm) < 1e-15)
  {
    const double_v_t e2bNorm = sqrt(sqr(b3) + sqr(b1));             // |e_2 x b|
    const double_v_t be2bNorm = sqrt(sqr(b1)*sqr(b2) + sqr(sqr(b1)*sqr(b3)) + sqr(b2) + sqr(b3));   // |b x (e_3 x b)|

    if (MathUtility::abs(e2bNorm) < 1e-15 || MathUtility::abs(be2bNorm) < 1e-15)
    {
      LOG(ERROR) << "Trying to transform a matrix to direction vector " << directionVector;
      return;
    }

    // use e_2

    Matrix<3,3,double_v_t> basis(
    {
      // normalized(b)    normalized(e_2 x b)   normalized(b x e_2 x b)
      b1/bNorm,            b3/e2bNorm,                       -b1*b2/be2bNorm,
      b2/bNorm,                     0,          (sqr(b1) + sqr(b3))/be2bNorm,
      b3/bNorm,           -b1/e2bNorm,                       -b2*b3/be2bNorm
    });

    Matrix<3,3,double_v_t> basis_T(
    {
      basis(0,0), basis(1,0), basis(2,0),
      basis(0,1), basis(1,1), basis(2,1),
      basis(0,2), basis(1,2), basis(2,2)
    });

    matrix = basis * matrix * basis_T;
  }
  else
  {
    // use e_3

    Matrix<3,3,double_v_t> basis(
    {
      // normalized(b)    normalized(e_3 x b)   normalized(b x e_3 x b)
      b1/bNorm,           -b2/e3bNorm,                       -b1*b3/be3bNorm,
      b2/bNorm,            b1/e3bNorm,                       -b2*b3/be3bNorm,
      b3/bNorm,                     0,          (sqr(b1) + sqr(b2))/be3bNorm
    });

    Matrix<3,3,double_v_t> basis_T(
    {
      basis(0,0), basis(1,0), basis(2,0),
      basis(0,1), basis(1,1), basis(2,1),
      basis(0,2), basis(1,2), basis(2,2)
    });

    matrix = basis * matrix * basis_T;
  }
}

template<typename T, std::size_t N>
bool containsNanOrInf(std::array<T,N> value)
{
  for (std::size_t i = 0; i < N; i++)
  {
    if (containsNanOrInf(value.at(i)))
    {
      return true;
    }
  }
  return false;
}

}  // namespace
