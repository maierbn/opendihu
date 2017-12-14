#pragma once

#include <array>
#include "control/types.h"

#include "basis_function/tensor_product_base.h"

namespace BasisFunction
{

// base class to compute global dof no
template<int D,typename BasisFunctionType>
class TensorProductDim : public TensorProductBase<D,BasisFunctionType>
{
};

// partial specialization for D=1
template<typename BasisFunctionType>
class TensorProductDim<1,BasisFunctionType> : public TensorProductBase<1,BasisFunctionType>
{
public:
  //! return the global dof number of element-local dof dofIndex of element elementNo, nElements is the total number of elements
  static int getDofNo(element_idx_t elementNo, int dofIndex, std::array<int,1> nElements);
  
  //! compute the jacobian matrix
  static std::array<Vec3,1> computeJacobian(std::array<Vec3,TensorProductBase<1,BasisFunctionType>::nNodesPerElement()> &nodes,
                                            std::array<double,1> xi);
  
};

// partial specialization for D=2
template<typename BasisFunctionType>
class TensorProductDim<2,BasisFunctionType> : public TensorProductBase<2,BasisFunctionType>
{
public:
 
  //! return the global dof number of element-local dof dofIndex of element elementNo, nElements is the total number of elements
  static int getDofNo(element_idx_t elementNo, int dofIndex, std::array<int,2> nElements);
  
  //! compute the jacobian matrix
  static std::array<Vec3,2> computeJacobian(std::array<Vec3,TensorProductBase<2,BasisFunctionType>::nNodesPerElement()> &nodes,
                                            std::array<double,2> xi);
  
};

// partial specialization for D=3
template<typename BasisFunctionType>
class TensorProductDim<3,BasisFunctionType> : public TensorProductBase<3,BasisFunctionType>
{
public:
 
  //! return the global dof number of element-local dof dofIndex of element elementNo, nElements is the total number of elements
  static int getDofNo(element_idx_t elementNo, int dofIndex, std::array<int,3> nElements);
  
  //! compute the jacobian matrix
  static std::array<Vec3,3> computeJacobian(std::array<Vec3,TensorProductBase<3,BasisFunctionType>::nNodesPerElement()> &nodes,
                                            std::array<double,3> xi);
  
};

}  // namespace

#include "basis_function/tensor_product_dim.tpp"