#include "basis_on_mesh/01_basis_on_mesh_function.h"

#include "utility/math_utility.h"
#include "easylogging++.h"

#include <cmath>
#include <array>

namespace BasisOnMesh
{

template<typename MeshType,typename BasisFunctionType>
int BasisOnMeshFunction<MeshType,BasisFunctionType>::
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
double BasisOnMeshFunction<MeshType,BasisFunctionType>::
phi(int dofIndex, std::array<double,MeshType::dim()> xi)
{
  double result = 1.0;
  for(int dimNo = 0; dimNo < MeshType::dim(); dimNo++)
  {
    int basisFunctionIndex1D = BasisOnMeshFunction<MeshType,BasisFunctionType>::getBasisFunctionIndex1D(dofIndex, dimNo);
    result *= BasisFunctionType::phi(basisFunctionIndex1D,xi[dimNo]);
  }
  return result;
}

template<typename MeshType,typename BasisFunctionType>
double BasisOnMeshFunction<MeshType,BasisFunctionType>::
dPhidxi(int dofIndex, int derivativeIdx, std::array<double,MeshType::dim()> xi)
{
  double result = 1.0;
  for(int dimNo = 0; dimNo < MeshType::dim(); dimNo++)
  {
    int basisFunctionIndex1D = BasisOnMeshFunction<MeshType,BasisFunctionType>::getBasisFunctionIndex1D(dofIndex, dimNo);
    if (dimNo == derivativeIdx)
      result *= BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[dimNo]);
    else
      result *= BasisFunctionType::phi(basisFunctionIndex1D, xi[dimNo]);
  }
  return result;
}

template<typename MeshType,typename BasisFunctionType>
std::array<double,MeshType::dim()> BasisOnMeshFunction<MeshType,BasisFunctionType>::
gradPhi(int dofIndex, std::array<double,MeshType::dim()> xi)
{
  std::array<double,MeshType::dim()> gradient;
  for(int gradientEntryNo = 0; gradientEntryNo < MeshType::dim(); gradientEntryNo++)
  {
    gradient[gradientEntryNo] = BasisOnMeshFunction<MeshType,BasisFunctionType>::dPhidxi(dofIndex, gradientEntryNo, xi);
  }
  return gradient;
}


};  // namespace