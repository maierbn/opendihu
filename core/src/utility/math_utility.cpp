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
Vc::double_v pow(Vc::double_v base, double exponent)
{
  return base.apply([exponent](double d)
  {
    return std::pow(d, exponent);
  });
}

template<>
double pow(double base, double exponent)
{
  return std::pow(base, exponent);
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

template<>
double normSquared<1>(const VecD<1,Vc::double_v> node)
{
  Vc::double_v vcResult = sqr(node[0]);
  double result = 0;
  for (int vcComponent = 0; vcComponent < Vc::double_v::size(); vcComponent++)
  {
    result += sqr(vcResult[vcComponent]);
  }
  return result;
}

template<>
double normSquared<2>(const VecD<2,Vc::double_v> node)
{
  Vc::double_v vcResult = sqr(node[0]) + sqr(node[1]);
  double result = 0;
  for (int vcComponent = 0; vcComponent < Vc::double_v::size(); vcComponent++)
  {
    result += sqr(vcResult[vcComponent]);
  }
  return result;
}

template<>
double normSquared<3>(const Vec3_v node)
{
  Vc::double_v vcResult = sqr(node[0]) + sqr(node[1]) + sqr(node[2]);
  double result = 0;
  for (int vcComponent = 0; vcComponent < Vc::double_v::size(); vcComponent++)
  {
    result += sqr(vcResult[vcComponent]);
  }
  return result;
}

//! arc cosine
template<>
double acos<double>(double value)
{
  return acos(value);
}

template<>
Vc::double_v acos<Vc::double_v>(Vc::double_v value)
{
  return value.apply([](double v){return std::acos(v);});
}

template<>
double abs<double>(double value)
{
  return fabs(value);
}

template<>
Vc::double_v abs<Vc::double_v>(Vc::double_v value)
{
  return Vc::abs(value);
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

//! transform a 3xD2 matrix to a DxD matrix by filling up with identity entries
template<>
std::array<Vec3_v,3> transformToDxD<3,3>(const std::array<Vec3_v,3> &matrix)
{
  return matrix;
}

template<>
std::array<Vec3_v,3> transformToDxD<3,2>(const std::array<Vec3_v,2> &matrix)
{
  return std::array<Vec3_v,3>({
    matrix[0],
    matrix[1],
    Vec3_v({Vc::double_v::Zero(), Vc::double_v::Zero(), Vc::double_v::One()})
  });
}

template<>
std::array<std::array<Vc::double_v,2>,2> transformToDxD<2,2>(const std::array<Vec3_v,2> &matrix)
{
  return std::array<std::array<Vc::double_v,2>,2>({
    std::array<Vc::double_v,2>({matrix[0][0], matrix[0][1]}),
    std::array<Vc::double_v,2>({matrix[1][0], matrix[1][1]})
  });
}

template<>
std::array<Vec3_v,3> transformToDxD<3,1>(const std::array<Vec3_v,1> &matrix)
{
  return std::array<Vec3_v,3>({
    matrix[0],
    Vec3_v({Vc::double_v::Zero(), Vc::double_v::One(), Vc::double_v::Zero()}),
    Vec3_v({Vc::double_v::Zero(), Vc::double_v::Zero(), Vc::double_v::One()})}
  );
}

template<>
std::array<std::array<Vc::double_v,1>,1> transformToDxD<1,1>(const std::array<Vec3_v,1> &matrix)
{
  return std::array<std::array<Vc::double_v,1>,1>({std::array<Vc::double_v,1>({matrix[0][0]})});
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

template<>
bool equals<1>(std::array<double,1> a, std::array<double,1> b, double tolerance)
{
  return fabs(a[0] - b[0]) <= tolerance;
}

template<>
bool equals<2>(std::array<double,2> a, std::array<double,2> b, double tolerance)
{
  return fabs(a[0] - b[0]) <= tolerance && fabs(a[1] - b[1]) <= tolerance;
}

template<>
bool equals<3>(std::array<double,3> a, std::array<double,3> b, double tolerance)
{
  return fabs(a[0] - b[0]) <= tolerance && fabs(a[1] - b[1]) <= tolerance && fabs(a[2] - b[2]) <= tolerance;
}

template<>
bool isFinite<double>(double value)
{
  return std::isfinite(value);
}

template<>
bool isFinite<Vc::double_v>(Vc::double_v value)
{
  return Vc::all_of(Vc::isfinite(value));
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
  // determine the two coordinates out of {x,y,z} that will be used
  int coordinate0 = 0;
  int coordinate1 = 1;

  // determine bounding box of element
  double xMin = geometryValues[0][0];
  double xMax = geometryValues[0][0];
  double yMin = geometryValues[0][1];
  double yMax = geometryValues[0][1];
  double zMin = geometryValues[0][2];
  double zMax = geometryValues[0][2];

  for (int i = 1; i < 4; i++)
  {
    xMin = std::min(xMin, geometryValues[i][0]);
    xMax = std::max(xMax, geometryValues[i][0]);
    yMin = std::min(yMin, geometryValues[i][1]);
    yMax = std::max(yMax, geometryValues[i][1]);
    zMin = std::min(zMin, geometryValues[i][2]);
    zMax = std::max(zMax, geometryValues[i][2]);
  }
  double extentX = fabs(xMax - xMin);
  double extentY = fabs(yMax - yMin);
  double extentZ = fabs(zMax - zMin);

  // the lowest extent will not be used as coordinate
  /*if (extentX < extentY && extentX < extentZ)
  {
    // use y and z
    coordinate0 = 1;
    coordinate1 = 2;
  }
  else if (extentY <= extentX && extentY < extentZ)
  {
    // use x and z
    coordinate0 = 0;
    coordinate1 = 2;
  }
  else*/
  {
    // use x and y
    coordinate0 = 0;
    coordinate1 = 1;
  }
    coordinate0 = 1;
    coordinate1 = 2;

  // derivation using sympy in script invert_mapping.py
  const double xp1 = point[coordinate0];
  const double xp2 = point[coordinate1];

  double x11 = geometryValues[0][coordinate0];
  double x12 = geometryValues[0][coordinate1];

  const double x21 = geometryValues[1][coordinate0];
  const double x22 = geometryValues[1][coordinate1];

  const double x31 = geometryValues[2][coordinate0];
  const double x32 = geometryValues[2][coordinate1];

  const double x41 = geometryValues[3][coordinate0];
  const double x42 = geometryValues[3][coordinate1];

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

  VLOG(1) << "quadrilateralGetPointCoordinates, extents: [" << extentX << "," << extentY << "," << extentZ
    << "], coordinates: [" << coordinate0 << "," << coordinate1 << "], element: [" << x11 << "," << x12
    << "], [" << x21 << "," << x22 << "], [" << x31 << "," << x32 << "], [" << x41 << "," << x42 << "], p: ["
    << xp1 << "," << xp2 << "], xi: [" << xi1 << "," << xi2 << "]";
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

    //VLOG(1) << " eigenvalue estimation, it " << i << ", v: " << v << ", value: " << normV << ", error: " << fabs(normVPrevious - normV);
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
