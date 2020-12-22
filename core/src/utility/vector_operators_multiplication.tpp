#include "utility/vector_operators.h"

#include <sstream>

#include "utility/petsc_utility.h"
#include "utility/string_utility.h"
#include "control/types.h"
#include "easylogging++.h"

//! vector*scalar multiplication
template<typename double_v_t, std::size_t nComponents>
std::array<double_v_t,nComponents> operator*(std::array<double_v_t,nComponents> vector, double lambda)
{
  std::array<double_v_t,nComponents> result;

  //#pragma omp simd
  for (int i = 0; i < nComponents; i++)
  {
    result[i] = lambda * vector[i];
  }
  return result;
}

//! vector*scalar multiplication
template<typename double_v_t, std::size_t nComponents>
std::array<double_v_t,nComponents> operator*(std::array<double_v_t,nComponents> vector, Vc::double_v lambda)
{
  std::array<double_v_t,nComponents> result;

  //#pragma omp simd
  for (int i = 0; i < nComponents; i++)
  {
    result[i] = lambda * vector[i];
  }
  return result;
}

//! vector/scalar division
template<typename double_v_t, std::size_t nComponents>
std::array<double_v_t,nComponents> operator/(std::array<double_v_t,nComponents> vector, Vc::double_v lambda)
{
  std::array<double_v_t,nComponents> result;

  //#pragma omp simd
  for (int i = 0; i < nComponents; i++)
  {
    result[i] = lambda / vector[i];
  }
  return result;
}

//! scalar*vector multiplication
template<typename double_v_t, std::size_t nComponents>
std::array<double_v_t,nComponents> operator*(Vc::double_v lambda, std::array<double_v_t,nComponents> vector)
{
  std::array<double_v_t,nComponents> result;

  //#pragma omp simd
  for (int i = 0; i < nComponents; i++)
  {
    result[i] = lambda * vector[i];
  }
  return result;
}

//! scalar*vector multiplication
template<typename double_v_t, std::size_t nComponents>
std::array<double_v_t,nComponents> operator*(double lambda, std::array<double_v_t,nComponents> vector)
{
  std::array<double_v_t,nComponents> result;

  //#pragma omp simd
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
template<std::size_t M, std::size_t N, typename double_v_t, typename double_v_t2>
std::array<double_v_t,M> operator*(const std::array<std::array<double_v_t,M>,N> &matrix, const std::array<double_v_t2,N> vector)
{
  std::array<double_v_t,M> result({0.0});

  // column index
  for (int j = 0; j < N; j++)
  {
    // row index
    for (int i = 0; i < M; i++)
    {
      result[i] += matrix[j][i] * vector[j];
    }
  }
  return result;
}
