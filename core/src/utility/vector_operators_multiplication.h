#pragma once

#include <Python.h>  // has to be the first included header
#include <array>
#include <map>
#include <set>
#include <vector>
#include <functional>
#include <petscmat.h>
#include <limits>
#include <vc_or_std_simd.h>  // this includes <Vc/Vc> or a Vc-emulating wrapper of <experimental/simd> if available

/** This file contains elemental operators for vectors, stored as `std::array<double,nComponents>`.
 *  The template typename double_v_t usually stands for double and Vc::double_v.
 */

//! scalar*vector multiplication
template<typename double_v_t, std::size_t nComponents>
std::array<double_v_t,nComponents> operator*(double lambda, std::array<double_v_t,nComponents> vector);

//! scalar*vector multiplication
template<typename double_v_t, std::size_t nComponents>
std::array<double_v_t,nComponents> operator*(Vc::double_v lambda, std::array<double_v_t,nComponents> vector);

//! vector*scalar multiplication
template<typename double_v_t, std::size_t nComponents>
std::array<double_v_t,nComponents> operator*(std::array<double_v_t,nComponents> vector, double lambda);

//! vector*scalar multiplication
template<typename double_v_t, std::size_t nComponents>
std::array<double_v_t,nComponents> operator*(std::array<double_v_t,nComponents> vector, Vc::double_v lambda);

//! vector/scalar division
template<typename double_v_t, std::size_t nComponents>
std::array<double_v_t,nComponents> operator/(std::array<double_v_t,nComponents> vector, Vc::double_v lambda);

//! vector*scalar multiplication
template<typename T>
std::vector<T> operator*(std::vector<T> vector, double lambda);

//! component-wise vector multiplication
template<std::size_t nComponents>
std::array<double,nComponents> operator*(std::array<double,nComponents> vector1, std::array<double,nComponents> vector2); // component-wise multiplication

//! vector multiplication, outer product
template<std::size_t nComponents1, std::size_t nComponents2>
std::array<std::array<double,nComponents1>,nComponents2> operator*(const std::array<double,nComponents2> vector1, const std::array<double,nComponents1> vector2);

//! matrix-vector multiplication, note that there is a matrix class with also matrix-vector multiplication. It stores matrices in row-major order, here column-major order is assumed
template<std::size_t M, std::size_t N, typename double_v_t, typename double_v_t2>
std::array<double_v_t,M> operator*(const std::array<std::array<double_v_t,M>,N> &matrix, const std::array<double_v_t2,N> vector);

#include "utility/vector_operators_multiplication.tpp"
