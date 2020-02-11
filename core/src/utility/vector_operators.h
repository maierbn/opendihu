#pragma once

#include <Python.h>  // has to be the first included header
#include <array>
#include <map>
#include <set>
#include <vector>
#include <functional>
#include <petscmat.h>

/** This file contains elemental operators for vectors, stored as `std::array<double,nComponents>`.
 */

//! arbitrary type difference
template<typename T, std::size_t nComponents>
std::array<T,nComponents> operator-(std::array<T,nComponents> vector1, std::array<T,nComponents> vector2);

//! arbitrary type addition
template<typename T, std::size_t nComponents>
std::array<T,nComponents> operator+(std::array<T,nComponents> vector1, std::array<T,nComponents> vector2);

//! vector unary minus
template<typename T, std::size_t nComponents>
std::array<T,nComponents> operator-(const std::array<T,nComponents> &vector1);

//! vector increment operation
template<typename T, std::size_t nComponents>
std::array<T,nComponents> &operator+=(std::array<T,nComponents> &vector1, std::array<T,nComponents> vector2);

//! vector multiply operation
template<std::size_t nComponents>
std::array<double,nComponents> &operator*=(std::array<double,nComponents> &vector1, double lambda);

//! vector division operation
template<std::size_t nComponents>
std::array<double,nComponents> &operator/=(std::array<double,nComponents> &vector1, double lambda);

//! scalar*vector multiplication
template<std::size_t nComponents>
std::array<double,nComponents> operator*(double lambda, std::array<double,nComponents> vector);

//! vector*scalar multiplication
template<std::size_t nComponents>
std::array<double,nComponents> operator*(std::array<double,nComponents> vector, double lambda);

//! vector*scalar multiplication
template<typename T>
std::vector<T> operator*(std::vector<T> vector, double lambda);

//! component-wise vector multiplication
template<std::size_t nComponents>
std::array<double,nComponents> operator*(std::array<double,nComponents> vector1, std::array<double,nComponents> vector2); // component-wise multiplication

//! vector multiplication, outer product
template<std::size_t nComponents1, std::size_t nComponents2>
std::array<std::array<double,nComponents1>,nComponents2> operator*(const std::array<double,nComponents2> vector1, const std::array<double,nComponents1> vector2);

//! component-wise division
template<typename T, std::size_t nComponents>
std::array<T,nComponents> operator/(std::array<T,nComponents> vector1, std::array<T,nComponents> vector2);

//! scalar division
template<typename T, std::size_t nComponents>
std::array<T,nComponents> operator/(std::array<T,nComponents> vector1, double value);

//! matrix-vector multiplication, note that there is a matrix class with also matrix-vector multiplication. It stores matrices in row-major order, here column-major order is assumed
template<std::size_t M, std::size_t N>
std::array<double,M> operator*(const std::array<std::array<double,M>,N> &matrix, const std::array<double,N> vector);

//! comparison operator with double value, true if any of the components fulfills the conditions " < value"
template<typename T, std::size_t N>
bool operator<(const std::array<T,N> &vector, double value);

//! output array content to stream
template<typename T, std::size_t N>
std::ostream &operator<<(std::ostream &stream, const std::array<T,N> &vector);

//! output array content to stream
template<std::size_t N>
std::ostream &operator<<(std::ostream &stream, const std::array<std::size_t,N> vector);

template<typename T>
std::ostream &operator<<(std::ostream &stream, std::reference_wrapper<T> value);

//! output arbitrary vector
template<typename T, typename A>
std::ostream &operator<<(std::ostream &stream, const std::vector<T,A> &vector);

template<>
std::ostream &operator<<(std::ostream &stream, const std::vector<double> &vector);

//! output contents of stringstream
//std::ostream &operator<<(std::ostream &stream, const std::stringstream &stringstream);

//! comparison operator for vectors of arbitrary type
template<typename T>
bool operator==(const std::vector<T> &vector1, const std::vector<T> &vector2);

//! output operator for pairs of arbitrary type
template<typename T1, typename T2>
std::ostream &operator<<(std::ostream &stream, const std::pair<T1,T2> &pair);

//! output operator for maps of arbitrary type
template<typename T1, typename T2>
std::ostream &operator<<(std::ostream &stream, const std::map<T1,T2> &map);

//! output operator for sets of arbitrary type
template<typename T>
std::ostream &operator<<(std::ostream &stream, const std::set<T> &set);

//! output operators for tuples or arbitrary type
template <size_t index, typename... T>
typename std::enable_if<(index >= sizeof...(T))>::type
  getString(std::ostream &stream, const std::tuple<T...> &tuple);

template <size_t index, typename... T>
typename std::enable_if<(index < sizeof...(T))>::type
  getString(std::ostream &stream, const std::tuple<T...> &tuple);

template <typename... T>
std::ostream &operator<<(std::ostream& stream, const std::tuple<T...> &tuple);

//! output operator for PETSc matrices
//std::ostream &operator<<(std::ostream &stream, const Mat &mat);

//! output operator for PETSc vectors
//std::ostream &operator<<(std::ostream &stream, const Vec &vec);

#include "utility/vector_operators.tpp"
