#pragma once

#include <Python.h>  // has to be the first included header
#include <array>

// define types
typedef signed int node_idx_t;       // type to hold value of node no
//typedef unsigned int node_idx_t;       // type to hold value of node no
typedef node_idx_t element_idx_t;   // type to hold value of element no
typedef node_idx_t dof_idx_t;       // type to hold value of dof no
typedef std::array<double, 2> Vec2;    // 2D vector
typedef std::array<double, 3> Vec3;    // 3D vector to store position of a node
