#pragma once

#include <Python.h>  // has to be the first included header
#include <array>
#include <petscmat.h>
#include <Vc/Vc>

//#define USE_VECTORIZED_FE_MATRIX_ASSEMBLY     // if the vectorized implementation of integrating the stiffness and mass matrices should be used
// This options reduces the runtime significantly, theoretically by a factor of 4 for matrix assembly. However, the compile time increases by a factor of circa 4 to 5.
// example: examples/solid_mechanics/dynamic_mooney_rivlin/muscle_with_fat runtime: 5:08.61 min -> 3:07.93 min, compile time: 2m28,546s+1m34,024s -> 2m16,100s+9m10,411s
// This option is set by the build system.

// define types for element, node and dof numbering, for the local numberings
typedef PetscInt node_no_t;               // type to hold value of node no or number of nodes
typedef node_no_t element_no_t;           // type to hold value of element no or number of elements
typedef node_no_t dof_no_t;               // type to hold value of dof no or number of dofs

typedef unsigned long long global_no_t;   // type for global numbers

// define vector types which are simply aliases for std::array<double,D>
template<int D, typename double_v_t = double>
using VecD = std::array<double_v_t,D>;    //< vector with D entries, actually the type of the 2nd template to std::array is of type std::size_t

using Vec2   = VecD<2>;    //< 2D vector
using Vec3   = VecD<3>;    //< 3D vector to store position of a node
using Vec3_v = VecD<3,Vc::double_v>;   //< SIMD vectorized types that use Vc::double_v instead of double and can be used to compute 4 values at once

// define tensor types
template<int D, typename double_v_t = double>
using Tensor2 = std::array<std::array<double_v_t,D>,D>;  //< two-point tensor of dimension DxD, storage is as vector of column vectors, i.e. column-major

template<int D, typename double_v_t = double>
using Tensor4 = std::array<std::array<std::array<std::array<double_v_t,D>,D>,D>,D>;  //< tensor of forth order of dimension DxD, storage is as vector of column vectors, i.e. column-major




// define data types,
// - either the vectorized version with Vc::double_v, which contains 4 double values. Then everything is computed for 4 consecutive elements at once
// - or use the normal, non-vectorized version. Then, normal double types are used.
#ifdef USE_VECTORIZED_FE_MATRIX_ASSEMBLY

// implementation with explicit SIMD vectorization
const int nVcComponents = Vc::double_v::size();   //< number of components that will be computed at once
using double_v_t = Vc::double_v;    //< data type for a vectorized value
using dof_no_v_t = Vc::int_v;       //< type for dof no indices, nVcComponents at once (actually more, but only the first nVcComponents are used)
using Vec3_v_t = Vec3_v;            //< instead of vectorized Vec3 use normal Vec3

template<int D>
using VecD_v_t = VecD<D,Vc::double_v>;         //< instead of vectorized VecD use a normal VecD

template<int D>
using Tensor2_v_t = Tensor2<D,Vc::double_v>;   //< instead of vectorized Tensor2 use a normal Tensor2

template<int D>
using Tensor4_v_t = Tensor4<D,Vc::double_v>;   //< instead of vectorized Tensor4 use a normal Tensor4

#else

// normal implementation without SIMD
const int nVcComponents = 1;   //< number of components that will be computed at once
using double_v_t = double;     //< data type for a normal value
using dof_no_v_t = dof_no_t;   //< type for dof no indices, nVcComponents at once (actually more, but only the first nVcComponents are used)
using Vec3_v_t = Vec3;            //< instead of vectorized Vec3 use normal Vec3

template<int D>
using VecD_v_t = VecD<D>;         //< instead of vectorized VecD use a normal VecD

template<int D>
using Tensor2_v_t = Tensor2<D>;   //< instead of vectorized Tensor2 use a normal Tensor2

template<int D>
using Tensor4_v_t = Tensor4<D>;   //< instead of vectorized Tensor4 use a normal Tensor4

#endif
