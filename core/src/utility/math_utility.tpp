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

template<int D>
double length(const VecD<D> node)
{
  return sqrt(sqr(node[0]) 
    + sqr(node[1])
    + sqr(node[2]));
}

}; // namespace
