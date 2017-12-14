#include "basis_function/tensor_product_base.h"

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
averageNDofsPerElement()
{
  return pow(BasisFunctionType::averageNDofsPerElement()/BasisFunctionType::nDofsPerNode(),D);
}

template<int D,typename BasisFunctionType>
constexpr int TensorProductBase<D,BasisFunctionType>::
nNodesPerElement()
{
  return pow(BasisFunctionType::nDofsPerBasis()/BasisFunctionType::nDofsPerNode(),D);
}


template<int D,typename BasisFunctionType>
int TensorProductBase<D,BasisFunctionType>::
getBasisFunctionIndex(int dofIndex, int dimNo)
{
  switch(dimNo)
  {
  case 0:
    return dofIndex % BasisFunctionType::nDofsPerBasis();
  case 1:
    return int((dofIndex % MathUtility::sqr(BasisFunctionType::nDofsPerBasis()))/ BasisFunctionType::nDofsPerBasis());
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
    int basisFunctionIndex = TensorProductBase<D,BasisFunctionType>::getBasisFunctionIndex(dofIndex, dimNo);
    result *= BasisFunctionType::phi(basisFunctionIndex,xi[dimNo]);
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
    gradient[gradientEntryNo] = 1.0;
    for(int dimNo = 0; dimNo < D; dimNo++)
    {
      int basisFunctionIndex = TensorProductBase<D,BasisFunctionType>::getBasisFunctionIndex(dofIndex, dimNo);
      if (dimNo == gradientEntryNo)
        gradient[gradientEntryNo] *= BasisFunctionType::dphi_dxi(basisFunctionIndex, xi[dimNo]);
      else
        gradient[gradientEntryNo] *= BasisFunctionType::phi(basisFunctionIndex, xi[dimNo]);
    }
  }
  return gradient;
}

};  // namespace