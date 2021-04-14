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

#ifdef HAVE_STDSIMD
  return std::experimental::pow(base, exponent);
#else
  return base.apply([exponent](double d)
  {
    return std::pow(d, exponent);
  });
#endif

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
    result += sqr((double)(vcResult[vcComponent]));
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
    result += sqr((double)(vcResult[vcComponent]));
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
    result += sqr((double)(vcResult[vcComponent]));
  }
  return result;
}

//! arc cosine
template<>
double acos<double>(double value)
{
  return std::acos(value);
}

template<>
Vc::double_v acos<Vc::double_v>(Vc::double_v value)
{
#ifdef HAVE_STDSIMD
  return std::experimental::acos(value);
#else
  return value.apply([](double v){return std::acos(v);});
#endif
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
    Vec3_v({Vc::double_v(Vc::Zero), Vc::double_v(Vc::Zero), Vc::double_v(Vc::One)})
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
    Vec3_v({Vc::double_v(Vc::Zero), Vc::double_v(Vc::One), Vc::double_v(Vc::Zero)}),
    Vec3_v({Vc::double_v(Vc::Zero), Vc::double_v(Vc::Zero), Vc::double_v(Vc::One)})}
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

void quadrilateralGetPointCoordinates(double xp1, double x11, double x21, double x31, double x41,
                                      double xp2, double x12, double x22, double x32, double x42,
                                      double &xi1, double &xi2)
{
  // Find the coordinates (xi1,xi2) such that xp = (1-xi1)*(1-xi2)*x1 + xi1*(1-xi2)*x2 + (1-xi1)*xi2*x3 + xi1*xi2*x4.
  // This formula was derived using sympy, in the file opendihu/doc/sympy/invert_mapping.py (section 2D linear quadrilateral).
  // There are two solutions to this problem, sometimes the correct solution is the first, sometimes it is the second.
  // Compute both and select the better one.

  // 1st solution
  const double xi1a = (1.0L/2.0L)*(2*x11*x32 - x11*x42 - 2*x12*x31 + x12*x41 - x21*x32 + x22*x31 + xp1*(x12 - x22 - x32 + x42) + xp2*(-x11 + x21 + x31 - x41) + std::sqrt(std::pow(x11, 2)*std::pow(x42, 2) - 2*std::pow(x11, 2)*x42*xp2 + std::pow(x11, 2)*std::pow(xp2, 2) - 2*x11*x12*x41*x42 + 2*x11*x12*x41*xp2 + 2*x11*x12*x42*xp1 - 2*x11*x12*xp1*xp2 - 2*x11*x21*x32*x42 + 2*x11*x21*x32*xp2 + 2*x11*x21*x42*xp2 - 2*x11*x21*std::pow(xp2, 2) - 2*x11*x22*x31*x42 + 2*x11*x22*x31*xp2 + 4*x11*x22*x32*x41 - 4*x11*x22*x32*xp1 - 4*x11*x22*x41*xp2 + 2*x11*x22*x42*xp1 + 2*x11*x22*xp1*xp2 + 2*x11*x31*x42*xp2 - 2*x11*x31*std::pow(xp2, 2) - 4*x11*x32*x41*xp2 + 2*x11*x32*x42*xp1 + 2*x11*x32*xp1*xp2 + 2*x11*x41*x42*xp2 + 2*x11*x41*std::pow(xp2, 2) - 2*x11*std::pow(x42, 2)*xp1 - 2*x11*x42*xp1*xp2 + std::pow(x12, 2)*std::pow(x41, 2) - 2*std::pow(x12, 2)*x41*xp1 + std::pow(x12, 2)*std::pow(xp1, 2) + 4*x12*x21*x31*x42 - 4*x12*x21*x31*xp2 - 2*x12*x21*x32*x41 + 2*x12*x21*x32*xp1 + 2*x12*x21*x41*xp2 - 4*x12*x21*x42*xp1 + 2*x12*x21*xp1*xp2 - 2*x12*x22*x31*x41 + 2*x12*x22*x31*xp1 + 2*x12*x22*x41*xp1 - 2*x12*x22*std::pow(xp1, 2) + 2*x12*x31*x41*xp2 - 4*x12*x31*x42*xp1 + 2*x12*x31*xp1*xp2 + 2*x12*x32*x41*xp1 - 2*x12*x32*std::pow(xp1, 2) - 2*x12*std::pow(x41, 2)*xp2 + 2*x12*x41*x42*xp1 - 2*x12*x41*xp1*xp2 + 2*x12*x42*std::pow(xp1, 2) + std::pow(x21, 2)*std::pow(x32, 2) - 2*std::pow(x21, 2)*x32*xp2 + std::pow(x21, 2)*std::pow(xp2, 2) - 2*x21*x22*x31*x32 + 2*x21*x22*x31*xp2 + 2*x21*x22*x32*xp1 - 2*x21*x22*xp1*xp2 + 2*x21*x31*x32*xp2 - 4*x21*x31*x42*xp2 + 2*x21*x31*std::pow(xp2, 2) - 2*x21*std::pow(x32, 2)*xp1 + 2*x21*x32*x41*xp2 + 2*x21*x32*x42*xp1 - 2*x21*x32*xp1*xp2 - 2*x21*x41*std::pow(xp2, 2) + 2*x21*x42*xp1*xp2 + std::pow(x22, 2)*std::pow(x31, 2) - 2*std::pow(x22, 2)*x31*xp1 + std::pow(x22, 2)*std::pow(xp1, 2) - 2*x22*std::pow(x31, 2)*xp2 + 2*x22*x31*x32*xp1 + 2*x22*x31*x41*xp2 + 2*x22*x31*x42*xp1 - 2*x22*x31*xp1*xp2 - 4*x22*x32*x41*xp1 + 2*x22*x32*std::pow(xp1, 2) + 2*x22*x41*xp1*xp2 - 2*x22*x42*std::pow(xp1, 2) + std::pow(x31, 2)*std::pow(xp2, 2) - 2*x31*x32*xp1*xp2 - 2*x31*x41*std::pow(xp2, 2) + 2*x31*x42*xp1*xp2 + std::pow(x32, 2)*std::pow(xp1, 2) + 2*x32*x41*xp1*xp2 - 2*x32*x42*std::pow(xp1, 2) + std::pow(x41, 2)*std::pow(xp2, 2) - 2*x41*x42*xp1*xp2 + std::pow(x42, 2)*std::pow(xp1, 2)))/(x11*x32 - x11*x42 - x12*x31 + x12*x41 - x21*x32 + x21*x42 + x22*x31 - x22*x41);
  const double xi2a = (x11*xi1a - x11 - x21*xi1a + xp1)/(x11*xi1a - x11 - x21*xi1a - x31*xi1a + x31 + x41*xi1a);

  // 2nd solution
  const double xi1b = (1.0L/2.0L)*(2*x11*x32 - x11*x42 - 2*x12*x31 + x12*x41 - x21*x32 + x22*x31 + xp1*(x12 - x22 - x32 + x42) + xp2*(-x11 + x21 + x31 - x41) - std::sqrt(std::pow(x11, 2)*std::pow(x42, 2) - 2*std::pow(x11, 2)*x42*xp2 + std::pow(x11, 2)*std::pow(xp2, 2) - 2*x11*x12*x41*x42 + 2*x11*x12*x41*xp2 + 2*x11*x12*x42*xp1 - 2*x11*x12*xp1*xp2 - 2*x11*x21*x32*x42 + 2*x11*x21*x32*xp2 + 2*x11*x21*x42*xp2 - 2*x11*x21*std::pow(xp2, 2) - 2*x11*x22*x31*x42 + 2*x11*x22*x31*xp2 + 4*x11*x22*x32*x41 - 4*x11*x22*x32*xp1 - 4*x11*x22*x41*xp2 + 2*x11*x22*x42*xp1 + 2*x11*x22*xp1*xp2 + 2*x11*x31*x42*xp2 - 2*x11*x31*std::pow(xp2, 2) - 4*x11*x32*x41*xp2 + 2*x11*x32*x42*xp1 + 2*x11*x32*xp1*xp2 + 2*x11*x41*x42*xp2 + 2*x11*x41*std::pow(xp2, 2) - 2*x11*std::pow(x42, 2)*xp1 - 2*x11*x42*xp1*xp2 + std::pow(x12, 2)*std::pow(x41, 2) - 2*std::pow(x12, 2)*x41*xp1 + std::pow(x12, 2)*std::pow(xp1, 2) + 4*x12*x21*x31*x42 - 4*x12*x21*x31*xp2 - 2*x12*x21*x32*x41 + 2*x12*x21*x32*xp1 + 2*x12*x21*x41*xp2 - 4*x12*x21*x42*xp1 + 2*x12*x21*xp1*xp2 - 2*x12*x22*x31*x41 + 2*x12*x22*x31*xp1 + 2*x12*x22*x41*xp1 - 2*x12*x22*std::pow(xp1, 2) + 2*x12*x31*x41*xp2 - 4*x12*x31*x42*xp1 + 2*x12*x31*xp1*xp2 + 2*x12*x32*x41*xp1 - 2*x12*x32*std::pow(xp1, 2) - 2*x12*std::pow(x41, 2)*xp2 + 2*x12*x41*x42*xp1 - 2*x12*x41*xp1*xp2 + 2*x12*x42*std::pow(xp1, 2) + std::pow(x21, 2)*std::pow(x32, 2) - 2*std::pow(x21, 2)*x32*xp2 + std::pow(x21, 2)*std::pow(xp2, 2) - 2*x21*x22*x31*x32 + 2*x21*x22*x31*xp2 + 2*x21*x22*x32*xp1 - 2*x21*x22*xp1*xp2 + 2*x21*x31*x32*xp2 - 4*x21*x31*x42*xp2 + 2*x21*x31*std::pow(xp2, 2) - 2*x21*std::pow(x32, 2)*xp1 + 2*x21*x32*x41*xp2 + 2*x21*x32*x42*xp1 - 2*x21*x32*xp1*xp2 - 2*x21*x41*std::pow(xp2, 2) + 2*x21*x42*xp1*xp2 + std::pow(x22, 2)*std::pow(x31, 2) - 2*std::pow(x22, 2)*x31*xp1 + std::pow(x22, 2)*std::pow(xp1, 2) - 2*x22*std::pow(x31, 2)*xp2 + 2*x22*x31*x32*xp1 + 2*x22*x31*x41*xp2 + 2*x22*x31*x42*xp1 - 2*x22*x31*xp1*xp2 - 4*x22*x32*x41*xp1 + 2*x22*x32*std::pow(xp1, 2) + 2*x22*x41*xp1*xp2 - 2*x22*x42*std::pow(xp1, 2) + std::pow(x31, 2)*std::pow(xp2, 2) - 2*x31*x32*xp1*xp2 - 2*x31*x41*std::pow(xp2, 2) + 2*x31*x42*xp1*xp2 + std::pow(x32, 2)*std::pow(xp1, 2) + 2*x32*x41*xp1*xp2 - 2*x32*x42*std::pow(xp1, 2) + std::pow(x41, 2)*std::pow(xp2, 2) - 2*x41*x42*xp1*xp2 + std::pow(x42, 2)*std::pow(xp1, 2)))/(x11*x32 - x11*x42 - x12*x31 + x12*x41 - x21*x32 + x21*x42 + x22*x31 - x22*x41);
  const double xi2b = (x11*xi1b - x11 - x21*xi1b + xp1)/(x11*xi1b - x11 - x21*xi1b - x31*xi1b + x31 + x41*xi1b);

  // choose the one that is closer to the interval [0,1]
  double scoreA = 0;
  if (xi1a < 0)
    scoreA = std::max(scoreA, -xi1a);
  else if (xi1a > 1)
    scoreA = std::max(scoreA, xi1a - 1);

  if (xi2a < 0)
    scoreA = std::max(scoreA, -xi2a);
  else if (xi2a > 1)
    scoreA = std::max(scoreA, xi2a - 1);

  double scoreB = 0;
  if (xi1b < 0)
    scoreB = std::max(scoreB, -xi1b);
  else if (xi1b > 1)
    scoreB = std::max(scoreB, xi1b - 1);

  if (xi2b < 0)
    scoreB = std::max(scoreB, -xi2b);
  else if (xi2b > 1)
    scoreB = std::max(scoreB, xi2b - 1);

  if (scoreA < scoreB)
  {
    xi1 = xi1a;
    xi2 = xi2a;
  }
  else
  {
    xi1 = xi1b;
    xi2 = xi2b;
  }

}


void quadrilateralGetPointCoordinates(const std::array<Vec3,4> geometryValues, const Vec3 point, Vec2 &xi)
{
  // derivation using sympy in script invert_mapping.py

  // project point onto plane of triangle
  // https://math.stackexchange.com/questions/544946/determine-if-projection-of-3d-point-onto-plane-is-within-a-triangle

  Vec3 u = geometryValues[1] - geometryValues[0];
  Vec3 v = geometryValues[2] - geometryValues[0];
  Vec3 w = point - geometryValues[0];
  Vec3 n = cross(u,v);
  double nSquared = normSquared<3>(n);
  double gamma = dot(cross(u,w), n) / nSquared;
  double beta = dot(cross(w,v), n) / nSquared;
  double alpha = 1 - gamma - beta;

  Vec3 pointProjected = alpha * geometryValues[0] + beta * geometryValues[1] + gamma * geometryValues[2];

  // We have 3D points but the problem whether the point is inside the quadrilateral is 2D, therefore an axis-aligned projection is performed
  // It is not clear which two coordinates should be used. Therefore we do all three choices and use the best result.

  // loop over different choices for the two coordinates to use
  std::array<int,3> coordinates0 = {0, 0, 1};
  std::array<int,3> coordinates1 = {1, 2, 2};
  std::array<int,3> coordinates2 = {2, 1, 0};
  double bestError = -1;

  for (int i = 0; i < coordinates0.size(); i++)
  {
    const int coordinate0 = coordinates0[i];
    const int coordinate1 = coordinates1[i];
    const int coordinate2 = coordinates2[i];

    // rename involved points
    const double xp1 = pointProjected[coordinate0];
    const double xp2 = pointProjected[coordinate1];

    const double x11 = geometryValues[0][coordinate0];
    const double x12 = geometryValues[0][coordinate1];

    const double x21 = geometryValues[1][coordinate0];
    const double x22 = geometryValues[1][coordinate1];

    const double x31 = geometryValues[2][coordinate0];
    const double x32 = geometryValues[2][coordinate1];

    const double x41 = geometryValues[3][coordinate0];
    const double x42 = geometryValues[3][coordinate1];

    // call formula
    double xi1, xi2;
    quadrilateralGetPointCoordinates(xp1, x11, x21, x31, x41,
                                     xp2, x12, x22, x32, x42, xi1, xi2);

    // validation, determine error
    const double test = (1-xi1)*(1-xi2)*geometryValues[0][coordinate2]
                        + xi1*(1-xi2)*geometryValues[1][coordinate2]
                        + (1-xi1)*xi2*geometryValues[2][coordinate2]
                        + xi1*xi2*geometryValues[3][coordinate2];
    double error = (test-pointProjected[coordinate2]) * (test-pointProjected[coordinate2]);

    if (bestError == -1 || error < bestError)
    {
      bestError = error;
      xi[0] = xi1;
      xi[1] = xi2;
    }
  }
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

bool containsNanOrInf(const double value)
{
  // check if the double value is a nan or a high value
  if (!std::isfinite(value) || fabs(value) > 1e+75 || value == std::numeric_limits<double>::max())
  {
    return true;
  }
  return false;
}

bool containsNanOrInf(const Vc::double_v value)
{
  for (int i = 0; i < Vc::double_v::size(); i++)
  {
    if (containsNanOrInf(value[i]))
    {
      return true;
    }
  }
  return false;
}

bool containsNanOrInf(const Vc::double_v value, dof_no_v_t elementNoLocalv)
{
  for (int i = 0; i < Vc::double_v::size(); i++)
  {
#ifdef USE_VECTORIZED_FE_MATRIX_ASSEMBLY
    // do not consider indices that are -1
    if (elementNoLocalv[i] == -1)
      continue;
#endif

    if (containsNanOrInf(value[i]))
    {
      return true;
    }
  }
  return false;
}

}  // namespace
