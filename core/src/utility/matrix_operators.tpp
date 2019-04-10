#include "utility/matrix_operators.h"


//! matrix difference
template<int nRows, int nColumns>
MathUtility::Matrix<nRows,nColumns> operator-(const MathUtility::Matrix<nRows,nColumns> matrix1, const MathUtility::Matrix<nRows,nColumns> matrix2)
{
  MathUtility::Matrix<nRows,nColumns> result;

  //#pragma omp simd
  for (int i = 0; i < nRows*nColumns; i++)
  {
    result[i] = matrix1[i] - matrix2[i];
  }
  return result;
}

//! matrix addition
template<int nRows, int nColumns>
MathUtility::Matrix<nRows,nColumns> operator+(const MathUtility::Matrix<nRows,nColumns> matrix1, const MathUtility::Matrix<nRows,nColumns> matrix2)
{
  MathUtility::Matrix<nRows,nColumns> result;

  //#pragma omp simd
  for (int i = 0; i < nRows*nColumns; i++)
  {
    result[i] = matrix1[i] + matrix2[i];
  }
  return result;
}

//! matrix increment operation
template<int nRows, int nColumns>
MathUtility::Matrix<nRows,nColumns> &operator+=(MathUtility::Matrix<nRows,nColumns> &matrix1, const MathUtility::Matrix<nRows,nColumns> matrix2)
{
  //#pragma omp simd
  for (int i = 0; i < nRows*nColumns; i++)
  {
    matrix1[i] += matrix2[i];
  }
  return matrix1;
}

//! scalar*matrix multiplication
template<int nRows, int nColumns>
MathUtility::Matrix<nRows,nColumns> operator*(double lambda, const MathUtility::Matrix<nRows,nColumns> matrix)
{
  MathUtility::Matrix<nRows,nColumns> result;

  //#pragma omp simd
  for (int i = 0; i < nRows*nColumns; i++)
  {
    result[i] = lambda * matrix[i];
  }
  return result;
}

//! matrix*scalar multiplication
template<int nRows, int nColumns>
MathUtility::Matrix<nRows,nColumns> operator*(MathUtility::Matrix<nRows,nColumns> matrix, double lambda)
{
  MathUtility::Matrix<nRows,nColumns> result;

  //#pragma omp simd
  for (int i = 0; i < nRows*nColumns; i++)
  {
    result[i] = lambda * matrix[i];
  }
  return result;
}

//! matrix-matrix multiplication
template<int nRows, int nColumns, int nColumns2>
MathUtility::Matrix<nRows,nColumns2> operator*(MathUtility::Matrix<nRows,nColumns> matrix1, MathUtility::Matrix<nColumns,nColumns2> matrix2)
{
  MathUtility::Matrix<nRows,nColumns2> result;

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


template<int nRows, int nColumns>
std::ostream &operator<<(std::ostream &stream, MathUtility::Matrix<nRows,nColumns> &matrix)
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
