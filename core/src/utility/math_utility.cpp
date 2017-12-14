#include "utility/math_utility.h"

#include <cmath>

double MathUtility::sqr(double v)
{
  return v*v;
}

int MathUtility::sqr(int v)
{
  return v*v;
}

double MathUtility::length(Vec3 node)
{
  return sqrt(MathUtility::sqr(node[0]) 
    + MathUtility::sqr(node[1])
    + MathUtility::sqr(node[2]));
  
}
double MathUtility::distance(Vec3 node1, Vec3 node2)
{
  return sqrt(MathUtility::sqr(node1[0]-node2[0]) 
    + MathUtility::sqr(node1[1]-node2[1])
    + MathUtility::sqr(node1[2]-node2[2]));
}

Vec3 operator-(Vec3 node1, Vec3 node2)
{
  Vec3 result;
  result[0] = node1[0]-node2[0];
  result[1] = node1[1]-node2[1];
  result[2] = node1[2]-node2[2];
  return result;
}

Vec3 operator+(Vec3 node1, Vec3 node2)
{
  Vec3 result;
  result[0] = node1[0]+node2[0];
  result[1] = node1[1]+node2[1];
  result[2] = node1[2]+node2[2];
  return result;
}

Vec3 operator*(double lambda, Vec3 node)
{
  Vec3 result;
  result[0] = lambda*node[0];
  result[1] = lambda*node[1];
  result[2] = lambda*node[2];
  return result;
}

std::array<double,9> MathUtility::computeTransformationMatrixAndDeterminant(Vec3& jacobianColumn0, Vec3& jacobianColumn1, 
                                                                            Vec3& jacobianColumn2, double &determinant)
{
  // rename input values
  double &m11 = jacobianColumn0[0];
  double &m21 = jacobianColumn0[1];
  double &m31 = jacobianColumn0[2];
  double &m12 = jacobianColumn1[0];
  double &m22 = jacobianColumn1[1];
  double &m32 = jacobianColumn1[2];
  double &m13 = jacobianColumn2[0];
  double &m23 = jacobianColumn2[1];
  double &m33 = jacobianColumn2[2];
  
  // the following code is copy&pasted from the output of ./invert_mapping.py 
  determinant = m11*m22*m33 - m11*m23*m32 - m12*m21*m33 + m12*m23*m31 + m13*m21*m32 - m13*m22*m31;
  double prefactor = 1./MathUtility::sqr(determinant);
  
  // 0 1 2
  // 3 4 5
  // 6 7 8
  
  // compute matrix entries of symmetric matrix
  std::array<double, 9> result;
  result[0] = prefactor * (MathUtility::sqr(m12*m23 - m13*m22) + MathUtility::sqr(m12*m33 - m13*m32) + MathUtility::sqr(m22*m33 - m23*m32));
  result[1] = prefactor * (-(m11*m23 - m13*m21)*(m12*m23 - m13*m22) - (m11*m33 - m13*m31)*(m12*m33 - m13*m32) - (m21*m33 - m23*m31)*(m22*m33 - m23*m32));
  result[2] = prefactor * ((m11*m22 - m12*m21)*(m12*m23 - m13*m22) + (m11*m32 - m12*m31)*(m12*m33 - m13*m32) + (m21*m32 - m22*m31)*(m22*m33 - m23*m32));
  result[3] = result[1];
  result[4] = prefactor * (MathUtility::sqr(m11*m23 - m13*m21) + MathUtility::sqr(m11*m33 - m13*m31) + MathUtility::sqr(m21*m33 - m23*m31));
  result[5] = prefactor * (-(m11*m22 - m12*m21)*(m11*m23 - m13*m21) - (m11*m32 - m12*m31)*(m11*m33 - m13*m31) - (m21*m32 - m22*m31)*(m21*m33 - m23*m31));
  result[6] = result[2];
  result[7] = result[5];
  result[8] = prefactor * (MathUtility::sqr(m11*m22 - m12*m21) + MathUtility::sqr(m11*m32 - m12*m31) + MathUtility::sqr(m21*m32 - m22*m31));
  
  return result;
}

double MathUtility::compute2DIntegrationFactor(Vec3& jacobianColumn0, Vec3& jacobianColumn1)
{
  // rename input values
  double &m11 = jacobianColumn0[0];
  double &m21 = jacobianColumn0[1];
  double &m31 = jacobianColumn0[2];
  double &m12 = jacobianColumn1[0];
  double &m22 = jacobianColumn1[1];
  double &m32 = jacobianColumn1[2];
  return sqrt(MathUtility::sqr(m11*m22) + MathUtility::sqr(m11*m32) - 2*m11*m12*m21*m22 - 2*m11*m12*m31*m32 
    + MathUtility::sqr(m12*m21) + MathUtility::sqr(m12*m31) + MathUtility::sqr(m21*m32) - 2*m21*m22*m31*m32 + MathUtility::sqr(m22*m31));
}

double MathUtility::applyTransformation(std::array<double,9> &transformationMatrix, Vec3 &vector1, Vec3 &vector2)
{
  double result;
  
  // rename input values
  double &v11 = vector1[0];
  double &v12 = vector1[1];
  double &v13 = vector1[2];
  double &v21 = vector2[0];
  double &v22 = vector2[1];
  double &v23 = vector2[2];
  // 0 1 2
  // 3 4 5
  // 6 7 8
  double &m11 = transformationMatrix[0];
  double &m12 = transformationMatrix[1];
  double &m13 = transformationMatrix[2];
  double &m22 = transformationMatrix[4];
  double &m23 = transformationMatrix[5];
  double &m33 = transformationMatrix[8];
  
  // compute result
  result = v21*(m11*v11 + m12*v12 + m13*v13) + v22*(m12*v11 + m22*v12 + m23*v13) + v23*(m13*v11 + m23*v12 + m33*v13);
  return result;
}

double MathUtility::applyTransformation(std::array<double,4> &transformationMatrix, Vec2 &vector1, Vec2 &vector2)
{
  double result;
  
  // rename input values
  double &v11 = vector1[0];
  double &v12 = vector1[1];
  double &v21 = vector2[0];
  double &v22 = vector2[1];
  // 0 1
  // 2 3
  double &m11 = transformationMatrix[0];
  double &m12 = transformationMatrix[1];
  double &m22 = transformationMatrix[3];
  
  // compute result
  result = v21*(m11*v11 + m12*v12) + v22*(m12*v11 + m22*v12);
  return result;
}

double MathUtility::computeDeterminant(Vec3& jacobianColumn0, Vec3& jacobianColumn1, Vec3& jacobianColumn2)
{
  // rename input values
  double &m11 = jacobianColumn0[0];
  double &m21 = jacobianColumn0[1];
  double &m31 = jacobianColumn0[2];
  double &m12 = jacobianColumn1[0];
  double &m22 = jacobianColumn1[1];
  double &m32 = jacobianColumn1[2];
  double &m13 = jacobianColumn2[0];
  double &m23 = jacobianColumn2[1];
  double &m33 = jacobianColumn2[2];
  
  return m11*m22*m33 - m11*m23*m32 - m12*m21*m33 + m12*m23*m31 + m13*m21*m32 - m13*m22*m31;
}
