#pragma once

#include <Python.h>  // has to be the first included header
#include <array>
#include "utility/matrix.h"

/** This file contains elemental operators for MathUtility::Matrix objects which are basically `std::array<double,nRows*nColumns>`.
 */

//! matrix difference
template<int nRows, int nColumns, typename double_v_t>
MathUtility::Matrix<nRows,nColumns,double_v_t> operator-(MathUtility::Matrix<nRows,nColumns,double_v_t> matrix1, MathUtility::Matrix<nRows,nColumns,double_v_t> matrix2);

//! matrix addition
template<int nRows, int nColumns, typename double_v_t>
MathUtility::Matrix<nRows,nColumns,double_v_t> operator+(MathUtility::Matrix<nRows,nColumns,double_v_t> matrix1, MathUtility::Matrix<nRows,nColumns,double_v_t> matrix2);

//! matrix increment operation
template<int nRows, int nColumns, typename double_v_t>
MathUtility::Matrix<nRows,nColumns,double_v_t> &operator+=(MathUtility::Matrix<nRows,nColumns,double_v_t> &matrix1, MathUtility::Matrix<nRows,nColumns,double_v_t> matrix2);

//! scalar*matrix multiplication
template<int nRows, int nColumns, typename double_v1_t, typename double_v2_t>
MathUtility::Matrix<nRows,nColumns,double_v2_t> operator*(double_v1_t lambda, MathUtility::Matrix<nRows,nColumns,double_v2_t> matrix);

//! matrix*scalar multiplication
template<int nRows, int nColumns, typename double_v1_t>
MathUtility::Matrix<nRows,nColumns,double_v1_t> operator*(MathUtility::Matrix<nRows,nColumns,double_v1_t> matrix, double_v1_t lambda);

//! matrix*scalar multiplication
template<int nRows, int nColumns, typename double_v1_t>
MathUtility::Matrix<nRows,nColumns,double_v1_t> operator*(MathUtility::Matrix<nRows,nColumns,double_v1_t> matrix, double lambda);

//! matrix-matrix multiplication
template<int nRows, int nColumns, int nColumns2, typename double_v_t>
MathUtility::Matrix<nRows,nColumns2,double_v_t> operator*(MathUtility::Matrix<nRows,nColumns,double_v_t> matrix1, MathUtility::Matrix<nColumns,nColumns2,double_v_t> matrix2);

//! output array content to stream
template<int nRows, int nColumns, typename double_v_t>
std::ostream &operator<<(std::ostream &stream, MathUtility::Matrix<nRows,nColumns,double_v_t> &matrix);

#include "utility/matrix_operators.tpp"
