#include "basis_function/tensor_product_dim.h"

#include <cmath>
#include <array>

namespace BasisFunction
{

// element-local dofIndex to global dofNo for 1D
template<typename BasisFunctionType>
int TensorProductDim<1,BasisFunctionType>::
getDofNo(element_idx_t elementNo, int dofIndex, std::array<int,1> nElements)
{
  // L linear  L quadratic  H quadratic
  // 0 1       0 1 2        0,1 2,3
  // averageNDofsPerElement: 
  // 1         2            2
  return BasisFunctionType::averageNDofsPerElement() * elementNo + dofIndex;
}

// element-local dofIndex to global dofNo for 2D
template<typename BasisFunctionType>
int TensorProductDim<2,BasisFunctionType>::
getDofNo(element_idx_t elementNo, int dofIndex, std::array<int,2> nElements)
{
  // L linear  quadratic  H quadratic
  //           6 7 8
  // 2 3       3 4 5      45 67
  // 0 1       0 1 2      01 23
  // nDofsPerBasis:
  // 2         3          4
  // averageNDofsPerElement:
  // 1         4          2
  int averageNDofsPerElement = BasisFunctionType::averageNDofsPerElement();
  int dofsPerRow = (averageNDofsPerElement * nElements[0] + 1);
  int elementX = int(elementNo % nElements[0]);
  int elementY = int(elementNo / nElements[0]);
  int localX = dofIndex % BasisFunctionType::nDofsPerBasis();
  int localY = int(dofIndex/BasisFunctionType::nDofsPerBasis());
  
  VLOG(2) << "  dof " << elementNo << ":" << dofIndex << ", element: ("<<elementX<<","<<elementY<<"), dofsPerRow="<<dofsPerRow<<", local: ("<<localX<<","<<localY<<")";
  
  return dofsPerRow * (elementY * averageNDofsPerElement + localY) 
    + averageNDofsPerElement * elementX + localX;
}

// element-local dofIndex to global dofNo for 3D
template<typename BasisFunctionType>
int TensorProductDim<3,BasisFunctionType>::
getDofNo(element_idx_t elementNo, int dofIndex, std::array<int,3> nElements)
{
  // L linear  quadratic  H quadratic
  //           6 7 8
  // 2 3       3 4 5      45 67
  // 0 1       0 1 2      01 23
  // nDofsPerBasis:
  // 2         3          4
  // averageNDofsPerElement:
  // 1         4          2
  int averageNDofsPerElement = BasisFunctionType::averageNDofsPerElement();
  int dofsPerRow = (averageNDofsPerElement * nElements[0] + 1);
  int dofsPerPlane = (averageNDofsPerElement * nElements[1] + 1) * dofsPerRow;
  
  int elementZ = int(elementNo / (nElements[0] * nElements[1]));
  int elementY = int((elementNo % (nElements[0] * nElements[1])) / nElements[0]);
  int elementX = elementNo % nElements[0];
  int localZ = int(dofIndex / MathUtility::sqr(BasisFunctionType::nDofsPerBasis()));
  int localY = int((dofIndex % MathUtility::sqr(BasisFunctionType::nDofsPerBasis())) / BasisFunctionType::nDofsPerBasis());
  int localX = dofIndex % BasisFunctionType::nDofsPerBasis();
  
  return dofsPerPlane * (elementZ * averageNDofsPerElement + localZ)
    + dofsPerRow * (elementY * averageNDofsPerElement + localY) 
    + averageNDofsPerElement * elementX + localX;
}

// Jacobian for 1D linear Lagrange basis
template<typename BasisFunctionType>
std::array<Vec3,1> TensorProductDim<1,BasisFunctionType>::
computeJacobian(std::array<Vec3,TensorProductBase<1,BasisFunctionType>::nNodesPerElement()> &node, std::array<double,1> xi)
{
  Vec3 jacobianColumn0 = (node[1]-node[0]);
  return std::array<Vec3,1>({jacobianColumn0});
}

// Jacobian for 2D linear Lagrange basis
template<typename BasisFunctionType>
std::array<Vec3,2> TensorProductDim<2,BasisFunctionType>::
computeJacobian(std::array<Vec3,TensorProductBase<2,BasisFunctionType>::nNodesPerElement()> &node, std::array<double,2> xi)
{
  double xi1 = xi[0];
  double xi2 = xi[1];
  Vec3 jacobianColumn0 = (1-xi2) * (node[1]-node[0]) + xi2 * (node[3]-node[2]);
  Vec3 jacobianColumn1 = (1-xi1) * (node[2]-node[0]) + xi1 * (node[3]-node[1]);
  return std::array<Vec3,2>({jacobianColumn0, jacobianColumn1});
}

// Jacobian for 3D linear Lagrange basis
template<typename BasisFunctionType>
std::array<Vec3,3> TensorProductDim<3,BasisFunctionType>::
computeJacobian(std::array<Vec3,TensorProductBase<3,BasisFunctionType>::nNodesPerElement()> &node, std::array<double,3> xi)
{
  double xi1 = xi[0];
  double xi2 = xi[1];
  double xi3 = xi[2];
  
  Vec3 jacobianColumn0 =
    (1-xi2) * (1-xi3) * (node[1]-node[0])
    + xi2 * (1-xi3) * (node[3]-node[2])
    + (1-xi2) * xi3 * (node[5]-node[4])
    + xi2 * xi3 * (node[7]-node[6]);

  Vec3 jacobianColumn1 =
    (1-xi1) * (1-xi3) * (node[2]-node[0])
    + xi1 * (1-xi3) * (node[3]-node[1])
    + (1-xi1) * xi3 * (node[6]-node[4])
    + xi1 * xi3 * (node[7]-node[5]);
  
  Vec3 jacobianColumn2 =
    (1-xi1) * (1-xi2) * (node[4]-node[0])
    + xi1 * (1-xi2) * (node[5]-node[1])
    + (1-xi1) * xi2 * (node[6]-node[2])
    + xi1 * xi2 * (node[7]-node[3]);
    
  return std::array<Vec3,3>({jacobianColumn0, jacobianColumn1, jacobianColumn2});
}


};  // namespace