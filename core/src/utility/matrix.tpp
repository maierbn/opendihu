#include "utility/matrix.h"

#pragma once

#include <Python.h>  // has to be the first included header
#include <array>
#include <algorithm>
#include <numeric>

namespace MathUtility 
{
  
template<int nRows, int nColumns>
Matrix<nRows,nColumns>::
Matrix(const std::array<double, nRows*nColumns> &rhs) : std::array<double, nRows*nColumns>(rhs)
{
}

template<int nRows, int nColumns>
Matrix<nRows,nColumns>::
Matrix(const std::array<double, nRows*nColumns> &&rhs) : std::array<double, nRows*nColumns>(rhs)
{
}

template<int nRows, int nColumns>
Matrix<nRows,nColumns>::
Matrix(double value) : std::array<double, nRows*nColumns>({value})
{
}

//! return a reference to the entry (rowIndex,columnIndex)
template<int nRows, int nColumns>
double &Matrix<nRows,nColumns>::
operator()(int rowIndex, int columnIndex)
{
  //return operator[](columnIndex*nRows + rowIndex);   // column-major
  return this->operator[](rowIndex*nColumns + columnIndex);   // row-major
}

//! fill a PETSc matrix with the own values
template<int nRows, int nColumns>
void Matrix<nRows,nColumns>::
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

template<int nRows, int nColumns>
std::array<double,nRows> Matrix<nRows,nColumns>::
operator*(const std::array<double,nColumns> &vector)
{
  std::array<double,nRows> result({0});
  for (int columnIndex = 0; columnIndex < nColumns; columnIndex++)
  {
#pragma simd
    for (int rowIndex = 0; rowIndex < nRows; rowIndex++)
    {
      result[rowIndex] += this->operator()(rowIndex, columnIndex) * vector[columnIndex];
    }
  }
  return result;
}

} // namespace
