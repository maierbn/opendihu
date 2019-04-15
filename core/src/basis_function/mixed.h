#pragma once

#include <Python.h>   // this has to be the first included header
#include <iostream>

namespace BasisFunction
{

/* class that holds the two basisfunctions for mixed formulation for solid mechanics
 */
template<typename LowOrderBasisFunctionType, typename HighOrderBasisFunctionType>
class Mixed
{
public:
  typedef LowOrderBasisFunctionType LowOrderBasisFunction;
  typedef HighOrderBasisFunctionType HighOrderBasisFunction;
  static const bool isMixed = true;  ///< if this is a mixed formulation comprising a higher order and a lower order basis function
};

template<typename BasisFunctionType>
using isNotMixed = std::enable_if_t<
  !BasisFunctionType::isMixed,
  BasisFunctionType
>;

}  // namespace
