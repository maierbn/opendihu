#pragma once

#include <Python.h>  // has to be the first included header
#include <array>
#include <petscmat.h>

// define types
//typedef signed int node_no_t;       // type to hold value of node no or number of nodes
//typedef std::size_t node_no_t;        // type to hold value of node no or number of nodes
// Warning: there is a problem when dof_no_t is incompatible to PetscInt, because this data type is used as index arrays in Vec and Mat functions

typedef PetscInt node_no_t;        // type to hold value of node no or number of nodes
typedef node_no_t element_no_t;       // type to hold value of element no or number of elements
typedef node_no_t dof_no_t;           // type to hold value of dof no or number of dofs

template<int D>
using VecD = std::array<double,D>;     // vector with D entries, actually the type of the 2nd template to std::array is of type std::size_t

typedef VecD<2> Vec2;    // 2D vector
typedef VecD<3> Vec3;    // 3D vector to store position of a node

template<int D>
using Tensor2 = std::array<std::array<double,D>,D>;  // two-point tensor of dimension DxD, storage is as vector of column vectors, i.e. column-major

template<int D>
using Tensor4 = std::array<std::array<std::array<std::array<double,D>,D>,D>,D>;  // tensor of forth order of dimension DxD, storage is as vector of column vectors, i.e. column-major

typedef unsigned long long global_no_t;   // type for global numbers
