#include "basis_function/tensor_product_dim.h"
//#include "basis_function/lagrange.h"

#include "easylogging++.h"

#include <cmath>
#include <array>

namespace BasisFunction
{

// element-local dofIndex to global dofNo for 1D
template<typename BasisFunctionType>
int TensorProductDim<1,BasisFunctionType>::
getDofNo(element_idx_t elementNo, int dofIndex, std::array<int,1> nElements)
{
  // L linear  L quadratic  H cubic
  // 0 1       0 1 2        0,1 2,3
  // averageNDofsPerElement: 
  // 1         2            2
  return TensorProductBase<1,BasisFunctionType>::averageNDofsPerElement() * elementNo + dofIndex;
}

// element-local dofIndex to global dofNo for 2D
template<typename BasisFunctionType>
int TensorProductDim<2,BasisFunctionType>::
getDofNo(element_idx_t elementNo, int dofIndex, std::array<int,2> nElements)
{
  // L linear  quadratic  H cubic
  //           6 7 8
  // 2 3       3 4 5      45 67
  // 0 1       0 1 2      01 23
  // nDofsPerBasis:
  // 2         3          4
  // averageNDofsPerElement:
  // 1         4          2
  int averageNDofsPerElement1D = TensorProductBase<1,BasisFunctionType>::averageNDofsPerElement();
  int dofsPerRow = (averageNDofsPerElement1D * nElements[0] + BasisFunctionType::nDofsPerNode());
  int elementX = int(elementNo % nElements[0]);
  int elementY = int(elementNo / nElements[0]);
  int localX = dofIndex % BasisFunctionType::nDofsPerBasis();
  int localY = int(dofIndex / BasisFunctionType::nDofsPerBasis());
  
  VLOG(2) << "  dof " << elementNo << ":" << dofIndex << ", element: ("<<elementX<<","<<elementY<<"), dofsPerRow="<<dofsPerRow<<", local: ("<<localX<<","<<localY<<")";
  
  return dofsPerRow * (elementY * averageNDofsPerElement1D + localY) 
    + averageNDofsPerElement1D * elementX + localX;
}

// element-local dofIndex to global dofNo for 3D
template<typename BasisFunctionType>
int TensorProductDim<3,BasisFunctionType>::
getDofNo(element_idx_t elementNo, int dofIndex, std::array<int,3> nElements)
{
  // L linear  quadratic  H cubic
  //           6 7 8
  // 2 3       3 4 5      45 67
  // 0 1       0 1 2      01 23
  // nDofsPerBasis:
  // 2         3          4
  // averageNDofsPerElement:
  // 1         4          2
  int averageNDofsPerElement1D = TensorProductBase<1,BasisFunctionType>::averageNDofsPerElement();
  int dofsPerRow = (averageNDofsPerElement1D * nElements[0] + BasisFunctionType::nDofsPerNode());
  int dofsPerPlane = (averageNDofsPerElement1D * nElements[1] + BasisFunctionType::nDofsPerNode()) * dofsPerRow;
  
  int elementZ = int(elementNo / (nElements[0] * nElements[1]));
  int elementY = int((elementNo % (nElements[0] * nElements[1])) / nElements[0]);
  int elementX = elementNo % nElements[0];
  int localZ = int(dofIndex / MathUtility::sqr(BasisFunctionType::nDofsPerBasis()));
  int localY = int((dofIndex % MathUtility::sqr(BasisFunctionType::nDofsPerBasis())) / BasisFunctionType::nDofsPerBasis());
  int localX = dofIndex % BasisFunctionType::nDofsPerBasis();
  
  return dofsPerPlane * (elementZ * averageNDofsPerElement1D + localZ)
    + dofsPerRow * (elementY * averageNDofsPerElement1D + localY) 
    + averageNDofsPerElement1D * elementX + localX;
}

// element-local nodeIndex to global nodeNo for 1D
template<typename BasisFunctionType>
int TensorProductDim<1,BasisFunctionType>::
getNodeNo(element_idx_t elementNo, int nodeIndex, std::array<int,1> nElements)
{
  // L linear  L quadratic  H cubic
  // 0 1       0 1 2        0 1
  // averageNDofsPerElement: 
  // 1         2            2
  // nDofsPerBasis:
  // 2         3            4
  // nDofsPerNode:
  // 1         1            2
  // nNodesPerElement:
  // 2         3            2
  return TensorProductBase<1,BasisFunctionType>::averageNNodesPerElement() * elementNo + nodeIndex;
}

// element-local nodeIndex to global nodeNo for 2D
template<typename BasisFunctionType>
int TensorProductDim<2,BasisFunctionType>::
getNodeNo(element_idx_t elementNo, int nodeIndex, std::array<int,2> nElements)
{
  // L linear  quadratic  H cubic
  //           6 7 8
  // 2 3       3 4 5      2 3
  // 0 1       0 1 2      0 1
  // nNodesPerElement:
  // 4         9          4
 
  int averageNNodesPerElement1D = TensorProductBase<1,BasisFunctionType>::averageNNodesPerElement();
  int nNodesPerElement1D = TensorProductBase<1,BasisFunctionType>::nNodesPerElement();
  int nodesPerRow = (averageNNodesPerElement1D * nElements[0] + 1);
  int elementX = int(elementNo % nElements[0]);
  int elementY = int(elementNo / nElements[0]);
  int localX = nodeIndex % nNodesPerElement1D;
  int localY = int(nodeIndex / nNodesPerElement1D);
  
  return nodesPerRow * (elementY * averageNNodesPerElement1D + localY) 
    + averageNNodesPerElement1D * elementX + localX;
}
// element-local nodeIndex to global nodeNo for 3D
template<typename BasisFunctionType>
int TensorProductDim<3,BasisFunctionType>::
getNodeNo(element_idx_t elementNo, int nodeIndex, std::array<int,3> nElements)
{
  int averageNNodesPerElement1D = TensorProductBase<1,BasisFunctionType>::averageNNodesPerElement();
  int nNodesPerElement1D = TensorProductBase<1,BasisFunctionType>::nNodesPerElement();
  int nodesPerRow = (averageNNodesPerElement1D * nElements[0] + 1);
  int nodesPerPlane = (averageNNodesPerElement1D * nElements[1] + 1) * nodesPerRow;
  
  int elementZ = int(elementNo / (nElements[0] * nElements[1]));
  int elementY = int((elementNo % (nElements[0] * nElements[1])) / nElements[0]);
  int elementX = elementNo % nElements[0];
  int localZ = int(nodeIndex / MathUtility::sqr(nNodesPerElement1D));
  int localY = int((nodeIndex % MathUtility::sqr(nNodesPerElement1D)) / nNodesPerElement1D);
  int localX = nodeIndex % nNodesPerElement1D;
  
  return nodesPerPlane * (elementZ * averageNNodesPerElement1D + localZ)
    + nodesPerRow * (elementY * averageNNodesPerElement1D + localY) 
    + averageNNodesPerElement1D * elementX + localX;
}

};  // namespace