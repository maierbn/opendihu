#pragma once

#include <Python.h>  // has to be the first included header
#include <array>

/** This file contains elemental operators for vectors, stored as `std::array<double,nComponents>`.
 */

//! vector difference
template<std::size_t nComponents>
std::array<double,nComponents> operator-(std::array<double,nComponents> vector1, std::array<double,nComponents> vector2);

//! vector addition
template<std::size_t nComponents>
std::array<double,nComponents> operator+(std::array<double,nComponents> vector1, std::array<double,nComponents> vector2);

//! vector increment operation
template<std::size_t nComponents>
std::array<double,nComponents> &operator+=(std::array<double,nComponents> &vector1, std::array<double,nComponents> vector2);

//! scalar*vector multiplication
template<std::size_t nComponents>
std::array<double,nComponents> operator*(double lambda, std::array<double,nComponents> vector);

//! vector*scalar multiplication
template<std::size_t nComponents>
std::array<double,nComponents> operator*(std::array<double,nComponents> vector, double lambda);

//! component-wise vector multiplication
template<std::size_t nComponents>
std::array<double,nComponents> operator*(std::array<double,nComponents> vector1, std::array<double,nComponents> vector2); // component-wise multiplication

//! output array content to stream
template<typename T, std::size_t N>
std::ostream &operator<<(std::ostream &stream, const std::array<T,N> &vector);

#if 0
//! output vector content to stream, already defined by SEMT library
template<typename T>
std::ostream &operator<<(std::ostream &stream, const std::vector<T> &vector);
#endif

//! comparison operator for vectors of arbitrary type
template<typename T>
bool operator==(const std::vector<T> &vector1, const std::vector<T> &vector2);

#include "utility/vector_operators.tpp"