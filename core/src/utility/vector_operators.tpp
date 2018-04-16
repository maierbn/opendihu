#include "utility/vector_operators.h"

#include "utility/petsc_utility.h"

//! vector difference
template<std::size_t nComponents>
std::array<double,nComponents> operator-(const std::array<double,nComponents> vector1, const std::array<double,nComponents> vector2)
{
  std::array<double,nComponents> result;
  
  #pragma simd
  for (int i = 0; i < nComponents; i++)
  {
    result[i] = vector1[i] - vector2[i];
  }
  return result;
}

//! vector addition
template<std::size_t nComponents>
std::array<double,nComponents> operator+(const std::array<double,nComponents> vector1, const std::array<double,nComponents> vector2)
{
  std::array<double,nComponents> result;
  
  #pragma simd
  for (int i = 0; i < nComponents; i++)
  {
    result[i] = vector1[i] + vector2[i];
  }
  return result;
}

//! vector increment operation
template<std::size_t nComponents>
std::array<double,nComponents> &operator+=(std::array<double,nComponents> &vector1, const std::array<double,nComponents> vector2)
{
  #pragma simd
  for (int i = 0; i < nComponents; i++)
  {
    vector1[i] += vector2[i];
  }
  return vector1;
}

//! scalar*vector multiplication
template<std::size_t nComponents>
std::array<double,nComponents> operator*(double lambda, const std::array<double,nComponents> vector)
{
  std::array<double,nComponents> result;
  
  #pragma simd
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
  
  #pragma simd
  for (int i = 0; i < nComponents; i++)
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
  
  #pragma simd
  for (int i = 0; i < nComponents; i++)
  {
    result[i] = vector1[i] * vector2[i];
  }
  return result;
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

#if 0
// operator already defined by SEMT
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
#endif


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