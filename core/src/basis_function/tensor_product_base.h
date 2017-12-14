#pragma once

#include <array>
#include "control/types.h"

namespace BasisFunction
{

template<int D,typename BasisFunctionType>
class TensorProductBase
{
public:
  //! number of degrees of freedom of this basis
  static constexpr int nDofsPerElement();
  
  //! if one assigns every dof to an element it is contained in, the number of degrees of freedom per element (not considering border elements)
  static constexpr int averageNDofsPerElement();
  
  //! number of nodes per element
  static constexpr int nNodesPerElement();
  
  //! evaluate the basis function corresponding to element-local dof dofIndex at xi, interval for xi is [0,1]
  static double phi(int dofIndex, std::array<double,D> xi);
  
  //! evaluate the first derivative of the basis function corresponding to element-local dof dofIndex at xi, interval for xi is [0,1]
  static std::array<double,D> gradPhi(int dofIndex, std::array<double,D> xi);

  //! return an array of all dof nos. of the element  
  static std::array<int,TensorProductBase<D,BasisFunctionType>::nDofsPerElement()> getElementDofs(element_idx_t elementNo, std::array<int,D> nElements);
  
private:
 
  //! given an element-local dofIndex and dimension No (0 <= dimNo < D), return the basis function index in that direction
  static int getBasisFunctionIndex(int dofIndex, int dimNo);
};


}  // namespace

#include "basis_function/tensor_product_base.tpp"