#include "utility/math_utility.h"

#include <cmath>

namespace MathUtility
{

template<int D>
double norm(const VecD<D> node)
{
  return length<D>(node);
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

}  // namespace
