#pragma once

#include <Python.h>  // has to be the first included header
#include <array>
//#include <Vc/Vc>

namespace MathUtility
{

/** Simple matrix helper class, the data is stored in row-major order
 */
template<int nRows, int nColumns>
struct Matrix :
  public std::array<double, nRows*nColumns>
  //public Vc::array<double, nRows*nColumns>
{
  //std::array<double, nRows*nColumns> data;   ///< the data member containing the entries in row-major order

  // -------- constructors -------------
  //! default empty constructor
  Matrix() = default;

  //! default copy constructor
  Matrix(const Matrix<nRows,nColumns> &) = default;

  //! copy constructor from array
  Matrix(const std::array<double, nRows*nColumns> &);

  //! default move constructor
  Matrix(Matrix<nRows,nColumns> &&) = default;

  //! move constructor from array
  Matrix(const std::array<double, nRows*nColumns> &&);

  //! default assignment operator
  Matrix<nRows,nColumns>& operator=(const Matrix<nRows,nColumns> &) = default;

  //! default move assignment operator
  Matrix<nRows,nColumns>& operator=(Matrix<nRows,nColumns> &&) = default;

  //! constructor to set all values at once (commented, because leads to error)
  //Matrix(double value);

  //! default destructor
  virtual ~Matrix() = default;


  // -------- other methods -------------
  //! return a reference to the entry (rowIndex,columnIndex)
  double &operator()(int rowIndex, int columnIndex);

  const double &operator()(int rowIndex, int columnIndex) const;

  //! fill a PETSc matrix with the own values
  void setPetscMatrix(Mat &mat);

  //! matrix*vector multiplication
  std::array<double,nRows> operator*(const std::array<double,nColumns> &vector);

};

} // namespace

#include "utility/matrix.tpp"
