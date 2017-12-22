#pragma once

#include <array>
#include "control/types.h"

#include "basis_function/tensor_product_base.h"
#include "basis_function/tensor_product_jacobian.h"

namespace BasisFunction
{

// base class to compute global dof no
template<int D,typename BasisFunctionType>
class TensorProductDim : public TensorProductBase<D,BasisFunctionType>
{
};

// partial specialization for D=1
template<typename BasisFunctionType>
class TensorProductDim<1,BasisFunctionType> : public TensorProductJacobian<1,BasisFunctionType>
{
public:
  //! return the global dof number of element-local dof dofIndex of element elementNo, nElements is the total number of elements
  static int getDofNo(element_idx_t elementNo, int dofIndex, std::array<int,1> nElements);
  
  //! return the global node number of element-local node nodeIndex of element elementNo, nElements is the total number of elements
  static int getNodeNo(element_idx_t elementNo, int nodeIndex, std::array<int,1> nElements);
};

// partial specialization for D=2
template<typename BasisFunctionType>
class TensorProductDim<2,BasisFunctionType> : public TensorProductJacobian<2,BasisFunctionType>
{
public:
 
  //! return the global dof number of element-local dof dofIndex of element elementNo, nElements is the total number of elements
  static int getDofNo(element_idx_t elementNo, int dofIndex, std::array<int,2> nElements);
  
  //! return the global node number of element-local node nodeIndex of element elementNo, nElements is the total number of elements
  static int getNodeNo(element_idx_t elementNo, int nodeIndex, std::array<int,2> nElements);
};

// partial specialization for D=3
template<typename BasisFunctionType>
class TensorProductDim<3,BasisFunctionType> : public TensorProductJacobian<3,BasisFunctionType>
{
public:
 
  //! return the global dof number of element-local dof dofIndex of element elementNo, nElements is the total number of elements
  static int getDofNo(element_idx_t elementNo, int dofIndex, std::array<int,3> nElements);
  
  //! return the global node number of element-local node nodeIndex of element elementNo, nElements is the total number of elements
  static int getNodeNo(element_idx_t elementNo, int nodeIndex, std::array<int,3> nElements);
};

}  // namespace

#include "basis_function/tensor_product_dim.tpp"