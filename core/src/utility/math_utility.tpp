#include "utility/math_utility.h"

template<typename T, unsigned long N>
std::ostream &operator<<(std::ostream &stream, const std::array<T,N> &node)
{
  stream << "(" << node[0];
  for (unsigned long i=1; i<N; i++)
    stream << "," << node[i];
  stream << ")";
  return stream;
}

template<typename T>
std::ostream &operator<<(std::ostream &stream, const std::vector<T> &values)
{
  if (values.empty())
  {
    stream << "()";
    return stream;
  }

  stream << "(" << values[0];
  for (unsigned long i=1; i<values.size(); i++)
    stream << "," << values[i];
  stream << ")";
  return stream;
}
