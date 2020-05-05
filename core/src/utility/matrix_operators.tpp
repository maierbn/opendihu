#include "utility/matrix_operators.h"


//! matrix difference
template<int nRows, int nColumns, typename double_v_t>
MathUtility::Matrix<nRows,nColumns,double_v_t> operator-(const MathUtility::Matrix<nRows,nColumns,double_v_t> matrix1, const MathUtility::Matrix<nRows,nColumns,double_v_t> matrix2)
{
  MathUtility::Matrix<nRows,nColumns,double_v_t> result;

  //#pragma omp simd
  for (int i = 0; i < nRows*nColumns; i++)
  {
    result[i] = matrix1[i] - matrix2[i];
  }
  return result;
}

//! matrix addition
template<int nRows, int nColumns, typename double_v_t>
MathUtility::Matrix<nRows,nColumns,double_v_t> operator+(const MathUtility::Matrix<nRows,nColumns,double_v_t> matrix1, const MathUtility::Matrix<nRows,nColumns,double_v_t> matrix2)
{
  MathUtility::Matrix<nRows,nColumns,double_v_t> result;

  //#pragma omp simd
  for (int i = 0; i < nRows*nColumns; i++)
  {
    result[i] = matrix1[i] + matrix2[i];
  }
  return result;
}

//! matrix increment operation
template<int nRows, int nColumns, typename double_v_t>
MathUtility::Matrix<nRows,nColumns,double_v_t> &operator+=(MathUtility::Matrix<nRows,nColumns,double_v_t> &matrix1, const MathUtility::Matrix<nRows,nColumns,double_v_t> matrix2)
{
  //#pragma omp simd
  for (int i = 0; i < nRows*nColumns; i++)
  {
    matrix1[i] += matrix2[i];
  }
  return matrix1;
}

//! scalar*matrix multiplication
template<int nRows, int nColumns, typename double_v1_t, typename double_v2_t>
MathUtility::Matrix<nRows,nColumns,double_v2_t> operator*(double_v1_t lambda, const MathUtility::Matrix<nRows,nColumns,double_v2_t> matrix)
{
  MathUtility::Matrix<nRows,nColumns,double_v2_t> result;

  //#pragma omp simd
  for (int i = 0; i < nRows*nColumns; i++)
  {
    result[i] = lambda * matrix[i];
  }
  return result;
}

//! matrix*scalar multiplication
template<int nRows, int nColumns, typename double_v1_t>
MathUtility::Matrix<nRows,nColumns,double_v1_t> operator*(MathUtility::Matrix<nRows,nColumns,double_v1_t> matrix, double_v1_t lambda)
{
  MathUtility::Matrix<nRows,nColumns,double_v1_t> result;

  //#pragma omp simd
  for (int i = 0; i < nRows*nColumns; i++)
  {
    result[i] = lambda * matrix[i];
  }
  return result;
}

// extra operator* when compiled with USE_VECTORIZED_FE_MATRIX_ASSEMBLY
#ifdef USE_VECTORIZED_FE_MATRIX_ASSEMBLY

//! matrix*scalar multiplication
template<int nRows, int nColumns, typename double_v1_t>
MathUtility::Matrix<nRows,nColumns,double_v1_t> operator*(MathUtility::Matrix<nRows,nColumns,double_v1_t> matrix, double lambda)
{
  MathUtility::Matrix<nRows,nColumns,double_v1_t> result;

  //#pragma omp simd
  for (int i = 0; i < nRows*nColumns; i++)
  {
    result[i] = lambda * matrix[i];
  }
  return result;
}

#endif

//! matrix-matrix multiplication
template<int nRows, int nColumns, int nColumns2, typename double_v_t>
MathUtility::Matrix<nRows,nColumns2,double_v_t> operator*(MathUtility::Matrix<nRows,nColumns,double_v_t> matrix1, MathUtility::Matrix<nColumns,nColumns2,double_v_t> matrix2)
{
  MathUtility::Matrix<nRows,nColumns2,double_v_t> result;

  //#pragma omp simd
  for (int i = 0; i < nRows; i++)
  {
    for (int j = 0; j < nColumns2; j++)
    {
      result[i*nColumns2 + j] = 0;

      for (int k = 0; k < nColumns; k++)
      {
        result[i*nColumns2 + j] += matrix1[i*nColumns + k] * matrix2[k*nColumns2 + j];
      }
    }
  }
  return result;
}


template<int nRows, int nColumns, typename double_v_t>
std::ostream &operator<<(std::ostream &stream, MathUtility::Matrix<nRows,nColumns,double_v_t> &matrix)
{
  if (matrix.empty())
  {
    stream << "()";
    return stream;
  }

  stream << "(" << matrix[0];
  for (unsigned long i = 1; i < matrix.size(); i++)
  {
    if ((i % nColumns) == 0)
    {
      stream << ";  ";
    }
    else
    {
      stream << ",";
    }
    stream << matrix[i];
  }
  stream << ")";
  return stream;
}
