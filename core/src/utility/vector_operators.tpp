#include "utility/vector_operators.h"

#include <sstream>

#include "utility/petsc_utility.h"
#include "easylogging++.h"

//! vector difference
template<typename T, std::size_t nComponents>
std::array<T,nComponents> operator-(const std::array<T,nComponents> vector1, const std::array<T,nComponents> vector2)
{
  std::array<T,nComponents> result;

  //#pragma omp simd
  for (int i = 0; i < nComponents; i++)
  {
    result[i] = vector1[i] - vector2[i];
  }
  return result;
}

//! vector unary minus
template<typename T, std::size_t nComponents>
std::array<T,nComponents> operator-(const std::array<T,nComponents> &vector1)
{
  std::array<T,nComponents> result;

  //#pragma omp simd
  for (int i = 0; i < nComponents; i++)
  {
    result[i] = -vector1[i];
  }
  return result;
}

//! vector addition
template<typename T, std::size_t nComponents>
std::array<T,nComponents> operator+(const std::array<T,nComponents> vector1, const std::array<T,nComponents> vector2)
{
  std::array<T,nComponents> result;

  //#pragma omp simd
  for (int i = 0; i < nComponents; i++)
  {
    result[i] = vector1[i] + vector2[i];
  }
  return result;
}

//! vector increment operation
template<typename T, std::size_t nComponents>
std::array<T,nComponents> &operator+=(std::array<T,nComponents> &vector1, const std::array<T,nComponents> vector2)
{
  //#pragma omp simd
  for (int i = 0; i < nComponents; i++)
  {
    vector1[i] += vector2[i];
  }
  return vector1;
}

//! vector multiply operation
template<std::size_t nComponents>
std::array<double,nComponents> &operator*=(std::array<double,nComponents> &vector1, double lambda)
{
  //#pragma omp simd
  for (int i = 0; i < nComponents; i++)
  {
    vector1[i] *= lambda;
  }
  return vector1;
}

//! vector division operation
template<std::size_t nComponents>
std::array<double,nComponents> &operator/=(std::array<double,nComponents> &vector1, double lambda)
{
  //#pragma omp simd
  for (int i = 0; i < nComponents; i++)
  {
    vector1[i] /= lambda;
  }
  return vector1;
}

//! scalar*vector multiplication
template<std::size_t nComponents>
std::array<double,nComponents> operator*(double lambda, const std::array<double,nComponents> vector)
{
  std::array<double,nComponents> result;

  //#pragma omp simd
  for (int i = 0; i < nComponents; i++)
  {
    result[i] = lambda * vector[i];
  }
  return result;
}

//! vector*scalar multiplication
template<std::size_t nComponents>
std::array<double,nComponents> operator*(std::array<double,nComponents> vector, double lambda)
{
  std::array<double,nComponents> result;

  //#pragma omp simd
  for (int i = 0; i < nComponents; i++)
  {
    result[i] = lambda * vector[i];
  }
  return result;
}

//! vector*scalar multiplication
template<typename T>
std::vector<T> operator*(std::vector<T> vector, double lambda)
{
  std::vector<T> result(vector.size());

  //#pragma omp simd
  for (int i = 0; i < vector.size(); i++)
  {
    result[i] = lambda * vector[i];
  }
  return result;
}

//! component-wise vector multiplication
template<std::size_t nComponents>
std::array<double,nComponents> operator*(const std::array<double,nComponents> vector1, const std::array<double,nComponents> vector2)
{
  std::array<double,nComponents> result;

  //#pragma omp simd
  for (int i = 0; i < nComponents; i++)
  {
    result[i] = vector1[i] * vector2[i];
  }
  return result;
}

//! vector multiplication, outer product
template<std::size_t nComponents1, std::size_t nComponents2>
std::array<std::array<double,nComponents1>,nComponents2> operator*(const std::array<double,nComponents2> vector1, const std::array<double,nComponents1> vector2)
{
  std::array<std::array<double,nComponents1>,nComponents2> result;

  //#pragma omp simd
  for (int i = 0; i < nComponents2; i++)
  {
    result[i] = vector1[i] * vector2;
  }
  return result;
}

//! matrix-vector multiplication
template<std::size_t M, std::size_t N>
std::array<double,M> operator*(const std::array<std::array<double,M>,N> &matrix, const std::array<double,N> vector)
{
  std::array<double,M> result({0.0});

  // column index
  for (int j = 0; j < N; j++)
  {
    // row index
    //#pragma omp simd
    for (int i = 0; i < M; i++)
    {
      result[i] += matrix[j][i] * vector[j];
    }
  }
  return result;
}

//! component-wise division
template<typename T, std::size_t nComponents>
std::array<T,nComponents> operator/(const std::array<T,nComponents> vector1, const std::array<T,nComponents> vector2)
{
  std::array<T,nComponents> result;

  //#pragma omp simd
  for (int i = 0; i < nComponents; i++)
  {
    result[i] = vector1[i] / vector2[i];
  }
  return result;
}

//! scalar division
template<typename T, std::size_t nComponents>
std::array<T,nComponents> operator/(const std::array<T,nComponents> vector1, const double value)
{
  std::array<T,nComponents> result;

  //#pragma omp simd
  for (int i = 0; i < nComponents; i++)
  {
    result[i] = vector1[i] / value;
  }
  return result;
}

template<typename T, std::size_t N>
bool operator<(const std::array<T,N> &vector, double value)
{
  for (int i = 0; i < N; i++)
  {
    if (vector[i] < value)
    {
      return true;
    }
  }

  return false;
}

//! output array content to stream
template<typename T, std::size_t N>
std::ostream &operator<<(std::ostream &stream, const std::array<T,N> &vector)
{
  stream << "(" << vector[0];
  for (std::size_t i = 1; i < N; i++)
    stream << "," << vector[i];
  stream << ")";
  return stream;
}

template<std::size_t N>
std::ostream &operator<<(std::ostream &stream, const std::array<std::size_t,N> vector)
{
  stream << "(" << vector[0];
  for (std::size_t i = 1; i < N; i++)
    stream << "," << vector[i];
  stream << ")";
  return stream;
}


template<typename T>
std::ostream &operator<<(std::ostream &stream, const std::vector<T> &values)
{
  if (values.empty())
  {
    stream << "[]";
    return stream;
  }

  stream << "[" << values[0];

  if (VLOG_IS_ON(1))
  {
    // with VLOG output all entries
    for (unsigned long i = 1; i < values.size(); i++)
      stream << "," << values[i];
  }
  else
  {
    // without VLOG only output the first 100 entries
    unsigned long i = 1;
    for (; i < std::min(100ul,values.size()); i++)
      stream << "," << values[i];
    if (i == 100 && i < values.size())
      stream << "... " << values.size() << " entries total, only showing the first 100 (set -vmodule=vector_operators*=1)";
  }

  stream << "]";
  return stream;
}

//! comparison operator for vectors of arbitrary type
template<typename T>
bool operator==(const std::vector<T> &vector1, const std::vector<T> &vector2)
{
  if (vector1.size() != vector2.size())
    return false;
  for (int i=0; i<vector1.size(); i++)
    if (vector1[i] != vector2[i])
      return false;
  return true;
}

//! output operator for pairs of arbitrary type
template<typename T1, typename T2>
std::ostream &operator<<(std::ostream &stream, const std::pair<T1,T2> &pair)
{
  stream << "(" << pair.first << "," << pair.second << ")";
  return stream;
}

//! output operator for maps of arbitrary type
template<typename T1, typename T2>
std::ostream &operator<<(std::ostream &stream, const std::map<T1,T2> &map)
{
  bool first = true;
  for (typename std::map<T1,T2>::const_iterator iter = map.cbegin(); iter != map.cend(); iter++)
  {
    if (!first)
      stream << ", ";
    stream << "\"" << iter->first << "\": " << iter->second;
    first = false;
  }
  return stream;
}

//! output operator for sets of arbitrary type
template<typename T>
std::ostream &operator<<(std::ostream &stream, const std::set<T> &set)
{
  stream << "{";
  bool first = true;
  for (typename std::set<T>::const_iterator iter = set.cbegin(); iter != set.cend(); iter++)
  {
    if (!first)
      stream << ", ";
    stream << (*iter);
    first = false;
  }
  stream << "}";
  return stream;
}

//! output operators for tuples or arbitrary type
template <size_t index, typename... T>
typename std::enable_if<(index >= sizeof...(T))>::type
  getString(std::ostream &stream, const std::tuple<T...> &tuple)
{}

template <size_t index, typename... T>
typename std::enable_if<(index < sizeof...(T))>::type
  getString(std::ostream &stream, const std::tuple<T...> &tuple)
{
  if (index != 0)
  {
    stream << ",";
  }
  stream << std::get<index>(tuple);

  getString<index+1>(stream, tuple);
}

template <typename... T>
std::ostream &operator<<(std::ostream& stream, const std::tuple<T...> &tuple)
{
  stream << "[";
  getString<0>(stream, tuple);
  stream << "]";

  return stream;
}

/*
std::ostream &operator<<(std::ostream &stream, const std::stringstream &stringstream)
{
  stream << "\"" << stringstream.str() << "\"";
  return stream;
}
*/
/*
std::ostream &operator<<(std::ostream &stream, const Mat &mat)
{
  int nRows, nColumns;
  MatGetSize(mat, &nRows, &nColumns);

  if (nRows*nColumns > 100)
  {
    stream << "Mat(" << nRows << "x" << nColumns << ")";
  }
  else
  {
    stream << PetscUtility::getStringMatrix(mat);
  }
  return stream;
}

std::ostream &operator<<(std::ostream &stream, const Vec &vec)
{
  int nEntries;
  VecGetSize(vec, &nEntries);

  if (nEntries > 100)
  {
    stream << "Vec(" << nEntries << ")";
  }
  else
  {
    stream << PetscUtility::getStringVector(vec);
  }
  return stream;
}*/
