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
};

}  // namespace