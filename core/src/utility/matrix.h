#pragma once

#include <Python.h>  // has to be the first included header
#include <array>
//#include <Vc/Vc>

namespace MathUtility
{

/** Simple matrix helper class, the data is stored in row-major order
 */
template<int nRows, int nColumns, typename double_v_t = double>
struct Matrix :
  public std::array<double_v_t, nRows*nColumns>
  //public Vc::array<double, nRows*nColumns>
{
  //std::array<double, nRows*nColumns> data;   ///< the data member containing the entries in row-major order

  // -------- constructors -------------
  //! default empty constructor
  Matrix() = default;

  //! default copy constructor
  Matrix(const Matrix<nRows,nColumns,double_v_t> &) = default;

  //! copy constructor from array
  Matrix(const std::array<double_v_t, nRows*nColumns> &);

  //! default move constructor
  Matrix(Matrix<nRows,nColumns,double_v_t> &&) = default;

  //! move constructor from array
  Matrix(const std::array<double_v_t, nRows*nColumns> &&);

  //! default assignment operator
  Matrix<nRows,nColumns,double_v_t>& operator=(const Matrix<nRows,nColumns,double_v_t> &) = default;

  //! default move assignment operator
  Matrix<nRows,nColumns,double_v_t>& operator=(Matrix<nRows,nColumns,double_v_t> &&) = default;

  //! constructor to set all values at once (commented, because leads to error)
  //Matrix(double_v_t value);

  //! default destructor
  virtual ~Matrix() = default;


  // -------- other methods -------------
  //! return a reference to the entry (rowIndex,columnIndex)
  double_v_t &operator()(int rowIndex, int columnIndex);

  const double_v_t &operator()(int rowIndex, int columnIndex) const;

  //! fill a PETSc matrix with the own values
  void setPetscMatrix(Mat &mat);

  //! matrix*vector multiplication
  template<typename double_v2_t>
  std::array<double_v_t,nRows> operator*(const std::array<double_v2_t,nColumns> &vector);

};

// type traits to find out if a type is a std::array<> or a Matrix<>
template<typename T>
struct isArrayOrMatrix : std::false_type {};

template<int nRows, int nColumns, typename double_v_t>
struct isArrayOrMatrix<Matrix<nRows,nColumns,double_v_t>> : std::true_type {};

template<typename double_v_t, std::size_t N>
struct isArrayOrMatrix<std::array<double_v_t,N>> : std::true_type {};

} // namespace

#include "utility/matrix.tpp"
