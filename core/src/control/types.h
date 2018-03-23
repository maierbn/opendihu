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
typedef std::array<double,2> Vec2;    // 2D vector, note that the type of the second template argument is std::size_t
typedef std::array<double,3> Vec3;    // 3D vector to store position of a node
