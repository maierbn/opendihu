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

double MathUtility::length(const Vec3 node)
{
  return sqrt(MathUtility::sqr(node[0]) 
    + MathUtility::sqr(node[1])
    + MathUtility::sqr(node[2]));
}

double MathUtility::norm(const Vec3 node)
{
  return MathUtility::length(node);
}
  
double MathUtility::distance(const Vec3 node1, const Vec3 node2)
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

Vec3 &operator+=(Vec3 &node1, Vec3 node2)
{
  node1[0] += node2[0];
  node1[1] += node2[1];
  node1[2] += node2[2];
  return node1;
}

std::ostream &operator<<(std::ostream &stream, std::array<double,1> node)
{
  stream << node[0];
  return stream;
}

template<typename T>
bool operator==(const std::vector<T> &vec1, const std::vector<T> &vec2)
{
  if (vec1.size() != vec2.size())
    return false;
  for (int i=0; i<vec1.size(); i++)
    if (vec1[i] != vec2[i])
      return false;
  return true;
}

std::array<double,9> MathUtility::computeTransformationMatrixAndDeterminant(const std::array<Vec3,3> &jacobian, double &determinant)
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

// 1D integration factor
template<>
double MathUtility::computeIntegrationFactor<1>(const std::array<Vec3,1> &jacobian)
{
  return MathUtility::length(jacobian[0]);
}

// 2D integration factor
template<>
double MathUtility::computeIntegrationFactor<2>(const std::array<Vec3,2> &jacobian)
{
  // rename input values
  const double m11 = jacobian[0][0];
  const double m21 = jacobian[0][1];
  const double m31 = jacobian[0][2];
  const double m12 = jacobian[1][0];
  const double m22 = jacobian[1][1];
  const double m32 = jacobian[1][2];
  return sqrt(MathUtility::sqr(m11*m22) + MathUtility::sqr(m11*m32) - 2*m11*m12*m21*m22 - 2*m11*m12*m31*m32 
    + MathUtility::sqr(m12*m21) + MathUtility::sqr(m12*m31) + MathUtility::sqr(m21*m32) - 2*m21*m22*m31*m32 + MathUtility::sqr(m22*m31));
}

// 3D integration factor
template<>
double MathUtility::computeIntegrationFactor<3>(const std::array<Vec3,3> &jacobian)
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

double MathUtility::applyTransformation(const std::array<double,9> &transformationMatrix, const Vec3 &vector1, const Vec3 &vector2)
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
  
  // compute result
  result = v21*(m11*v11 + m12*v12 + m13*v13) + v22*(m12*v11 + m22*v12 + m23*v13) + v23*(m13*v11 + m23*v12 + m33*v13);
  return result;
}

double MathUtility::applyTransformation(const std::array<double,4> &transformationMatrix, const Vec2 &vector1, const Vec2 &vector2)
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
  
  // compute result
  result = v21*(m11*v11 + m12*v12) + v22*(m12*v11 + m22*v12);
  return result;
}

double MathUtility::computeDeterminant(const std::array<Vec3,3> &jacobian)
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

bool MathUtility::isSubsequenceOf(std::vector<int> a, std::vector<int> b, int &subsequenceAStartPos)
{
  if (b.empty())
    return true;
  
  // find the matching entry in vector a
  bool matchFound = false;
  int aIndex=0;
  for (; aIndex<a.size(); aIndex++)
  {
    if (a[aIndex] == b[0])
    {
      subsequenceStartPos = aIndex;
      matchFound = true;
      break;
    }
  }

  if (!matchFound)
    return false;
  
  for (int bIndex=1; bIndex<b.size(); bIndex++)
  {
    aIndex++;
    if (aIndex >= a.size())
      return false;
    
    if (a[aIndex] != b[bIndex])
      return false;
  }
  return true;
}