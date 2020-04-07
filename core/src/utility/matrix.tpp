#include "utility/matrix.h"

#pragma once

#include <Python.h>  // has to be the first included header
#include <array>
#include <algorithm>
#include <numeric>

namespace MathUtility
{

template<int nRows, int nColumns, typename double_v_t>
Matrix<nRows,nColumns,double_v_t>::
Matrix(const std::array<double_v_t, nRows*nColumns> &rhs) : std::array<double_v_t, nRows*nColumns>(rhs)
{
}

template<int nRows, int nColumns, typename double_v_t>
Matrix<nRows,nColumns,double_v_t>::
Matrix(const std::array<double_v_t, nRows*nColumns> &&rhs) : std::array<double_v_t, nRows*nColumns>(rhs)
{
}
/*
template<int nRows, int nColumns, typename double_v_t>
Matrix<nRows,nColumns,double_v_t>::
Matrix(double_v_t value) : std::array<double_v_t, nRows*nColumns>({value})
{
}*/

//! return a reference to the entry (rowIndex,columnIndex)
template<int nRows, int nColumns, typename double_v_t>
double_v_t &Matrix<nRows,nColumns,double_v_t>::
operator()(int rowIndex, int columnIndex)
{
  assert(rowIndex >= 0);
  assert(rowIndex < nRows);
  assert(columnIndex >= 0);
  assert(columnIndex < nColumns);

  //return operator[](columnIndex*nRows + rowIndex);   // column-major
  return this->operator[](rowIndex*nColumns + columnIndex);   // row-major
}

//! return a const reference to the entry (rowIndex,columnIndex)
template<int nRows, int nColumns, typename double_v_t>
const double_v_t &Matrix<nRows,nColumns,double_v_t>::
operator()(int rowIndex, int columnIndex) const
{
  assert(rowIndex >= 0);
  assert(rowIndex < nRows);
  assert(columnIndex >= 0);
  assert(columnIndex < nColumns);

  //return operator[](columnIndex*nRows + rowIndex);   // column-major
  return this->operator[](rowIndex*nColumns + columnIndex);   // row-major
}

//! fill a PETSc matrix with the own values
template<int nRows, int nColumns, typename double_v_t>
void Matrix<nRows,nColumns,double_v_t>::
setPetscMatrix(Mat &mat)
{
  // prepare index arrays for MatSetValues
  std::array<int,nRows> rowIndices;
  std::iota(rowIndices.begin(), rowIndices.end(), 0);

  std::array<int,nColumns> columnIndices;
  std::iota(columnIndices.begin(), columnIndices.end(), 0);

  // assign values to PETSc data structure
  MatSetValues(mat, nRows, rowIndices.data, nColumns, columnIndices.data(), this->data(), INSERT_VALUES);
}

template<int nRows, int nColumns, typename double_v_t>
std::array<double_v_t,nRows> Matrix<nRows,nColumns,double_v_t>::
operator*(const std::array<double_v_t,nColumns> &vector)
{
  std::array<double_v_t,nRows> result({0});
  for (int columnIndex = 0; columnIndex < nColumns; columnIndex++)
  {
//#pragma omp simd
    for (int rowIndex = 0; rowIndex < nRows; rowIndex++)
    {
      result[rowIndex] += this->operator()(rowIndex, columnIndex) * vector[columnIndex];
    }
  }
  return result;
}

} // namespace
