#include "function_space/01_function_space_function.h"

#include "utility/math_utility.h"
#include "easylogging++.h"

#include <cmath>
#include <array>

namespace FunctionSpace
{

template<typename MeshType,typename BasisFunctionType>
int FunctionSpaceFunction<MeshType,BasisFunctionType>::
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

template<typename MeshType,typename BasisFunctionType>
double FunctionSpaceFunction<MeshType,BasisFunctionType>::
phi(int dofIndex, std::array<double,MeshType::dim()> xi)
{
  double result = 1.0;
  for (int dimNo = 0; dimNo < MeshType::dim(); dimNo++)
  {
    int basisFunctionIndex1D = FunctionSpaceFunction<MeshType,BasisFunctionType>::getBasisFunctionIndex1D(dofIndex, dimNo);
    result *= BasisFunctionType::phi(basisFunctionIndex1D,xi[dimNo]);
  }
  return result;
}

template<typename MeshType,typename BasisFunctionType>
double FunctionSpaceFunction<MeshType,BasisFunctionType>::
dphi_dxi(int dofIndex, int derivativeIdx, std::array<double,MeshType::dim()> xi)
{
  double result = 1.0;
  for (int dimNo = 0; dimNo < MeshType::dim(); dimNo++)
  {
    int basisFunctionIndex1D = FunctionSpaceFunction<MeshType,BasisFunctionType>::getBasisFunctionIndex1D(dofIndex, dimNo);
    if (dimNo == derivativeIdx)
      result *= BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[dimNo]);
    else
      result *= BasisFunctionType::phi(basisFunctionIndex1D, xi[dimNo]);
  }
  return result;
}

template<typename MeshType,typename BasisFunctionType>
std::array<double,MeshType::dim()> FunctionSpaceFunction<MeshType,BasisFunctionType>::
gradPhi(int dofIndex, std::array<double,MeshType::dim()> xi)
{
  std::array<double,MeshType::dim()> gradient;
  for (int gradientEntryNo = 0; gradientEntryNo < MeshType::dim(); gradientEntryNo++)
  {
    gradient[gradientEntryNo] = FunctionSpaceFunction<MeshType,BasisFunctionType>::dphi_dxi(dofIndex, gradientEntryNo, xi);
  }
  return gradient;
}

// for complete polynomials
template<typename MeshType, int order>
std::array<double,MeshType::dim()> FunctionSpaceFunction<MeshType, BasisFunction::CompletePolynomialOfDimensionAndOrder<MeshType::dim(),order>>::
gradPhi(int dofIndex, std::array<double,MeshType::dim()> xi)
{
  std::array<double,MeshType::dim()> gradient;
  for (int gradientEntryNo = 0; gradientEntryNo < MeshType::dim(); gradientEntryNo++)
  {
    gradient[gradientEntryNo] = BasisFunction::CompletePolynomialOfDimensionAndOrder<MeshType::dim(),order>::dphi_dxi(dofIndex, gradientEntryNo, xi);
  }
  return gradient;
}

};  // namespace
