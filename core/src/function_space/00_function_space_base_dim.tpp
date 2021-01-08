#include "function_space/00_function_space_base_dim.h"

#include "utility/math_utility.h"
#include "easylogging++.h"
#include "basis_function/hermite.h"
#include "basis_function/lagrange.h"

#include <cmath>
#include <array>

namespace FunctionSpace
{

// lagrange and hermite polynomials
template<int D,typename BasisFunctionType>
constexpr int FunctionSpaceBaseDim<D,BasisFunctionType>::
nDofsPerElement()
{
  return MathUtility::powConst<dof_no_t,int>(BasisFunctionType::nDofsPerBasis(),D);
}

template<int D,typename BasisFunctionType>
constexpr int FunctionSpaceBaseDim<D,BasisFunctionType>::
nNodesPerElement()
{
  return MathUtility::powConst<dof_no_t,int>(BasisFunctionType::nDofsPerBasis()/BasisFunctionType::nDofsPerNode(),D);
}

template<int D,typename BasisFunctionType>
constexpr int FunctionSpaceBaseDim<D,BasisFunctionType>::
nDofsPerNode()
{
  return nDofsPerElement() / nNodesPerElement();
}

template<int D,typename BasisFunctionType>
constexpr int FunctionSpaceBaseDim<D,BasisFunctionType>::
averageNNodesPerElement()
{
  return MathUtility::powConst<dof_no_t,int>(BasisFunctionType::nDofsPerBasis()/BasisFunctionType::nDofsPerNode()-1,D);
}

template<int D,typename BasisFunctionType>
constexpr int FunctionSpaceBaseDim<D,BasisFunctionType>::
averageNDofsPerElement()
{
 // nNodesPerBasis = nDofsPerBasis / nDofsPerNode
 // averageNNodesPerBasis = nNodesPerBasis - 1
 // averageNDofsPerElement = averageNNodesPerBasis * nDofsPerNode
  return FunctionSpaceBaseDim<D,BasisFunctionType>::averageNNodesPerElement() * BasisFunctionType::nDofsPerNode();
}

} // namespace
