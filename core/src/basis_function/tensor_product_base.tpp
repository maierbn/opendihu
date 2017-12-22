#include "basis_function/tensor_product_base.h"

#include "utility/math_utility.h"
#include "easylogging++.h"

#include <cmath>
#include <array>

namespace BasisFunction
{
  
template<int D,typename BasisFunctionType>
constexpr int TensorProductBase<D,BasisFunctionType>::
nDofsPerElement()
{
  return pow(BasisFunctionType::nDofsPerBasis(),D);
}

template<int D,typename BasisFunctionType>
constexpr int TensorProductBase<D,BasisFunctionType>::
nNodesPerElement()
{
  return pow(BasisFunctionType::nDofsPerBasis()/BasisFunctionType::nDofsPerNode(),D);
}

template<int D,typename BasisFunctionType>
constexpr int TensorProductBase<D,BasisFunctionType>::
averageNNodesPerElement()
{
  return pow(BasisFunctionType::nDofsPerBasis()/BasisFunctionType::nDofsPerNode()-1,D);
}

template<int D,typename BasisFunctionType>
constexpr int TensorProductBase<D,BasisFunctionType>::
averageNDofsPerElement()
{
 // nNodesPerBasis = nDofsPerBasis / nDofsPerNode
 // averageNNodesPerBasis = nNodesPerBasis - 1
 // averageNDofsPerElement = averageNNodesPerBasis * nDofsPerNode
  return TensorProductBase<D,BasisFunctionType>::averageNNodesPerElement() * BasisFunctionType::nDofsPerNode();
}

template<int D,typename BasisFunctionType>
int TensorProductBase<D,BasisFunctionType>::
getBasisFunctionIndex1D(int dofIndex, int dimNo)
{
  switch(dimNo)
  {
  case 0:
    return dofIndex % BasisFunctionType::nDofsPerBasis();
  case 1:
    return int((dofIndex % MathUtility::sqr(BasisFunctionType::nDofsPerBasis())) / BasisFunctionType::nDofsPerBasis());
  case 2:
    return int(dofIndex / MathUtility::sqr(BasisFunctionType::nDofsPerBasis()));
  default:
    return 0;
  }
}

template<int D,typename BasisFunctionType>
double TensorProductBase<D,BasisFunctionType>::
phi(int dofIndex, std::array<double,D> xi)
{
  double result = 1.0;
  for(int dimNo = 0; dimNo < D; dimNo++)
  {
    int basisFunctionIndex1D = TensorProductBase<D,BasisFunctionType>::getBasisFunctionIndex1D(dofIndex, dimNo);
    result *= BasisFunctionType::phi(basisFunctionIndex1D,xi[dimNo]);
  }
  return result;
}

template<int D,typename BasisFunctionType>
double TensorProductBase<D,BasisFunctionType>::
dPhidxi(int dofIndex, int derivativeIdx, std::array<double,D> xi)
{
  double result = 1.0;
  for(int dimNo = 0; dimNo < D; dimNo++)
  {
    int basisFunctionIndex1D = TensorProductBase<D,BasisFunctionType>::getBasisFunctionIndex1D(dofIndex, dimNo);
    if (dimNo == derivativeIdx)
      result *= BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[dimNo]);
    else
      result *= BasisFunctionType::phi(basisFunctionIndex1D, xi[dimNo]);
  }
  return result;
}

template<int D,typename BasisFunctionType>
std::array<double,D> TensorProductBase<D,BasisFunctionType>::
gradPhi(int dofIndex, std::array<double,D> xi)
{
  std::array<double,D> gradient;
  for(int gradientEntryNo = 0; gradientEntryNo < D; gradientEntryNo++)
  {
    gradient[gradientEntryNo] = TensorProductBase<D,BasisFunctionType>::dPhidxi(dofIndex, gradientEntryNo, xi);
  }
  return gradient;
}


};  // namespace