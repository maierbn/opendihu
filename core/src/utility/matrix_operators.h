#pragma once

#include <Python.h>  // has to be the first included header
#include <array>
#include "utility/matrix.h"

/** This file contains elemental operators for MathUtility::Matrix objects which are basically `std::array<double,nRows*nColumns>`.
 */

//! matrix difference
template<int nRows, int nColumns>
MathUtility::Matrix<nRows,nColumns> operator-(MathUtility::Matrix<nRows,nColumns> matrix1, MathUtility::Matrix<nRows,nColumns> matrix2);

//! matrix addition
template<int nRows, int nColumns>
MathUtility::Matrix<nRows,nColumns> operator+(MathUtility::Matrix<nRows,nColumns> matrix1, MathUtility::Matrix<nRows,nColumns> matrix2);

//! matrix increment operation
template<int nRows, int nColumns>
MathUtility::Matrix<nRows,nColumns> &operator+=(MathUtility::Matrix<nRows,nColumns> &matrix1, MathUtility::Matrix<nRows,nColumns> matrix2);

//! scalar*matrix multiplication
template<int nRows, int nColumns>
MathUtility::Matrix<nRows,nColumns> operator*(double lambda, MathUtility::Matrix<nRows,nColumns> matrix);

//! matrix*scalar multiplication
template<int nRows, int nColumns>
MathUtility::Matrix<nRows,nColumns> operator*(MathUtility::Matrix<nRows,nColumns> matrix, double lambda);

//! matrix-matrix multiplication
template<int nRows, int nColumns, int nColumns2>
MathUtility::Matrix<nRows,nColumns2> operator*(MathUtility::Matrix<nRows,nColumns> matrix1, MathUtility::Matrix<nColumns,nColumns2> matrix2);

//! output array content to stream
template<int nRows, int nColumns>
std::ostream &operator<<(std::ostream &stream, MathUtility::Matrix<nRows,nColumns> &matrix);

#include "utility/matrix_operators.tpp"
