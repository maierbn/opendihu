#include "utility/matrix.h"

#pragma once

#include <Python.h>  // has to be the first included header
#include <array>
#include <algorithm>

namespace MathUtility 
{
  
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

} // namespace
