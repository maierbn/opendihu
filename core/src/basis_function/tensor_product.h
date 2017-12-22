#pragma once

#include <array>
#include "control/types.h"
#include "basis_function/tensor_product_dim.h"

namespace BasisFunction
{

template<int D,typename BasisFunctionType>
class TensorProduct : public TensorProductDim<D,BasisFunctionType>
{
public:
  //! return an array of all dof nos. of the element  
  static std::array<int,TensorProductBase<D,BasisFunctionType>::nDofsPerElement()> getElementDofs(element_idx_t elementNo, std::array<int,D> nElements);

  //! return an array of all node nos. of the element  
  static std::array<int,TensorProductBase<D,BasisFunctionType>::nNodesPerElement()> getElementNodes(element_idx_t elementNo, std::array<int,D> nElements);

  //! return an array of the gradients of all nodal basis functions, evaluated at xi  
  static std::array<std::array<double,D>,TensorProductBase<D,BasisFunctionType>::nDofsPerElement()> getGradPhi(std::array<double,D> xi);
private:
};


}  // namespace

#include "basis_function/tensor_product.tpp"