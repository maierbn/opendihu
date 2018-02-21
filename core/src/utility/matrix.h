#pragma once

#include <Python.h>  // has to be the first included header
#include <array>

namespace MathUtility 
{

/** Simple matrix helper class
 */
template<int nRows, int nColumns>
struct Matrix :
  public std::array<double, nRows*nColumns>
{
  //std::array<double, nRows*nColumns> data;   ///< the data member containing the entries in row-major order
  
  //! return a reference to the entry (rowIndex,columnIndex)
  double &operator()(int rowIndex, int columnIndex);
  
  //! fill a PETSc matrix with the own values
  void setPetscMatrix(Mat &mat);
};

} // namespace

#include "utility/matrix.tpp"
