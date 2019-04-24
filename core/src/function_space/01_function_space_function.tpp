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
  const int nDofsPerNode = FunctionSpaceBaseDim<MeshType::dim(),BasisFunctionType>::nDofsPerNode();
  const int nDofsPerNode1D = BasisFunctionType::nDofsPerNode();
  int nNodesPerElement1D = FunctionSpaceBaseDim<1,BasisFunctionType>::nNodesPerElement();

  VLOG(2) << "D: " << MeshType::dim() << ", nDofsPerNode: " << nDofsPerNode << ", nNodesPerElement1D: " << nNodesPerElement1D;

  const int nodeNo = (int)(dofIndex/nDofsPerNode);
  const int nodalDofIndex = dofIndex % nDofsPerNode;    // nodal dof index, counting the D dimensional dofs

  int nodeIndexDimension = 0;
  int basisFunctionIndexPerNodeDimension = 0;   // this is always 0 for non-Hermite (because nodalDofIndex is 0, because nDofsPerNode is 1)

  switch(dimNo)
  {
  case 0:
    basisFunctionIndexPerNodeDimension = nodalDofIndex % nDofsPerNode1D;
    nodeIndexDimension = nodeNo % nNodesPerElement1D;
    break;
  case 1:
    basisFunctionIndexPerNodeDimension = int((nodalDofIndex % MathUtility::sqr(nDofsPerNode1D)) / nDofsPerNode1D);
    nodeIndexDimension = int((nodeNo % MathUtility::sqr(nNodesPerElement1D)) / nNodesPerElement1D);
    break;
  case 2:
    basisFunctionIndexPerNodeDimension = int(nodalDofIndex / MathUtility::sqr(nDofsPerNode1D));
    nodeIndexDimension = int(nodeNo / MathUtility::sqr(nNodesPerElement1D));
    break;
  default:
    LOG(ERROR) << "dimNo is not 0,1 or 2: " << dimNo;
    return 0;
  }

  VLOG(1) << "getBasisFunctionIndex1D(" << dofIndex << ", dimNo " << dimNo << "), nodalDofIndex: " << nodalDofIndex << ", nodeIndexDimension: " << nodeIndexDimension << ", basisFunctionIndexPerNodeDimension: " << basisFunctionIndexPerNodeDimension;

  return nodeIndexDimension * nDofsPerNode1D + basisFunctionIndexPerNodeDimension;
}

template<typename MeshType,typename BasisFunctionType>
double FunctionSpaceFunction<MeshType,BasisFunctionType>::
phi(int dofIndex, std::array<double,MeshType::dim()> xi)
{
  //VLOG(3) << "  -> phi(dofIndex " << dofIndex << ", xi " << xi << "), dim: " << MeshType::dim();
  double result = 1.0;
  for (int dimNo = 0; dimNo < MeshType::dim(); dimNo++)
  {
    int basisFunctionIndex1D = FunctionSpaceFunction<MeshType,BasisFunctionType>::getBasisFunctionIndex1D(dofIndex, dimNo);
    //VLOG(3) << "       (dim " << dimNo << ", xi=" << xi[dimNo] << ", basisFunctionIndex1D: " << basisFunctionIndex1D << ", phi: "
    //  << BasisFunctionType::phi(basisFunctionIndex1D,xi[dimNo]);
    result *= BasisFunctionType::phi(basisFunctionIndex1D,xi[dimNo]);
  }
  return result;
}

template<typename MeshType,typename BasisFunctionType>
double FunctionSpaceFunction<MeshType,BasisFunctionType>::
dphi_dxi(int dofIndex, int derivativeIdx, std::array<double,MeshType::dim()> xi)
{
  VLOG(1) << "   dphi" << dofIndex << "/dxi" << derivativeIdx << "(" << xi << ")";
  double result = 1.0;
  for (int dimNo = 0; dimNo < MeshType::dim(); dimNo++)
  {
    int basisFunctionIndex1D = FunctionSpaceFunction<MeshType,BasisFunctionType>::getBasisFunctionIndex1D(dofIndex, dimNo);
    if (dimNo == derivativeIdx)
    {
      result *= BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[dimNo]);
      VLOG(1) << "     dimNo " << dimNo << ", basisFunctionIndex1D " << basisFunctionIndex1D << "   result *= dphi_dxi [" << BasisFunctionType::dphi_dxi(basisFunctionIndex1D, xi[dimNo]) << "]";
    }
    else
    {
      result *= BasisFunctionType::phi(basisFunctionIndex1D, xi[dimNo]);
      VLOG(1) << "     dimNo " << dimNo << ", basisFunctionIndex1D " << basisFunctionIndex1D << "   result *= phi [" << BasisFunctionType::phi(basisFunctionIndex1D, xi[dimNo]) << "]";
    }
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

} // namespace
