#include "basis_on_mesh/04_basis_on_mesh_numbers_structured.h"

#include <Python.h>  // has to be the first included header
#include <cmath>
#include <array>

#include "easylogging++.h"

namespace BasisOnMesh
{

// These numberings run over the whole locally stored information, i.e. there is no distinguising between ghost and interior nodes here.
 
// element-local dofIndex to local dofNo for 1D
template<typename MeshType,typename BasisFunctionType>
dof_no_t BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getDofNo(element_no_t elementNo, int dofIndex) const
{
  // L linear  L quadratic  H cubic
  // 0 1       0 1 2        0,1 2,3
  // averageNDofsPerElement:
  // 1         2            2
  return BasisOnMeshFunction<MeshType,BasisFunctionType>::averageNDofsPerElement() * elementNo + dofIndex;
}

//! get all dofs of a specific node for 1D
template<typename MeshType,typename BasisFunctionType>
void BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const
{
  dofGlobalNos.reserve(dofGlobalNos.size() + BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerNode());
  for (int i=0; i<BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerNode(); i++)
  {
    dofGlobalNos.push_back(BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerNode() * nodeGlobalNo + i);
  }
}

//! get all dofs of a specific node for 1D
template<typename MeshType,typename BasisFunctionType>
void BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getNodeDofs(node_no_t nodeGlobalNo, std::array<dof_no_t,BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerNode()> &dofGlobalNos) const
{
  for (int i=0; i<BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerNode(); i++)
  {
    dofGlobalNos[i] = BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerNode() * nodeGlobalNo + i;
  }
}

//! get the dof no of the specified dof at the node for 1D
template<typename MeshType,typename BasisFunctionType>
dof_no_t BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getNodeDofNo(node_no_t nodeGlobalNo, int dofIndex) const
{
  return BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerNode() * nodeGlobalNo + dofIndex;
}

// element-local dofIndex to local dofNo for 2D
template<typename MeshType,typename BasisFunctionType>
dof_no_t BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
getDofNo(element_no_t elementNo, int dofIndex) const
{
  // L linear  quadratic  H cubic
  //           6 7 8
  // 2 3       3 4 5      45 67
  // 0 1       0 1 2      01 23
  // nDofsPerBasis:
  // 2         3          4
  // averageNDofsPerElement:
  // 1         4          2

  const std::array<element_no_t, MeshType::dim()> &nElements = this->nElementsPerCoordinateDirection_;
  int averageNDofsPerElement1D = BasisOnMeshBaseDim<1,BasisFunctionType>::averageNDofsPerElement();
  dof_no_t dofsPerRow = (averageNDofsPerElement1D * nElements[0] + BasisFunctionType::nDofsPerNode());
  element_no_t elementX = element_no_t(elementNo % nElements[0]);
  element_no_t elementY = element_no_t(elementNo / nElements[0]);
  dof_no_t localX = dofIndex % BasisFunctionType::nDofsPerBasis();
  dof_no_t localY = dof_no_t(dofIndex / BasisFunctionType::nDofsPerBasis());

  VLOG(4) << "  dof " << elementNo << ":" << dofIndex << ", element: ("<<elementX<<","<<elementY<<"), dofsPerRow="<<dofsPerRow<<", local: ("<<localX<<","<<localY<<")";

  return dofsPerRow * (elementY * averageNDofsPerElement1D + localY)
    + averageNDofsPerElement1D * elementX + localX;
}

//! get all dofs of a specific node for 2D
template<typename MeshType,typename BasisFunctionType>
void BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const
{
  dofGlobalNos.reserve(dofGlobalNos.size() + BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerNode());
  for (int i=0; i<BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerNode(); i++)
  {
    dofGlobalNos.push_back(BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerNode() * nodeGlobalNo + i);
  }
}

//! get all dofs of a specific node for 2D
template<typename MeshType,typename BasisFunctionType>
void BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
getNodeDofs(node_no_t nodeGlobalNo, std::array<dof_no_t,BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerNode()> &dofGlobalNos) const
{
  for (int i=0; i<BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerNode(); i++)
  {
    dofGlobalNos[i] = BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerNode() * nodeGlobalNo + i;
  }
}

//! get the dof no of the specified dof at the node for 2D
template<typename MeshType,typename BasisFunctionType>
dof_no_t BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
getNodeDofNo(node_no_t nodeGlobalNo, int dofIndex) const
{
  return BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerNode() * nodeGlobalNo + dofIndex;
}

// element-local dofIndex to local dofNo for 3D
template<typename MeshType,typename BasisFunctionType>
dof_no_t BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getDofNo(element_no_t elementNo, int dofIndex) const
{
  // L linear  quadratic  H cubic
  //           6 7 8
  // 2 3       3 4 5      45 67
  // 0 1       0 1 2      01 23
  // nDofsPerBasis:
  // 2         3          4
  // averageNDofsPerElement:
  // 1         4          2

  const std::array<element_no_t, MeshType::dim()> &nElements = this->nElementsPerCoordinateDirection_;
  int averageNDofsPerElement1D = BasisOnMeshBaseDim<1,BasisFunctionType>::averageNDofsPerElement();
  dof_no_t dofsPerRow = (averageNDofsPerElement1D * nElements[0] + BasisFunctionType::nDofsPerNode());
  dof_no_t dofsPerPlane = (averageNDofsPerElement1D * nElements[1] + BasisFunctionType::nDofsPerNode()) * dofsPerRow;

  element_no_t elementZ = element_no_t(elementNo / (nElements[0] * nElements[1]));
  element_no_t elementY = element_no_t((elementNo % (nElements[0] * nElements[1])) / nElements[0]);
  element_no_t elementX = elementNo % nElements[0];
  dof_no_t localZ = dof_no_t(dofIndex / MathUtility::sqr(BasisFunctionType::nDofsPerBasis()));
  dof_no_t localY = dof_no_t((dofIndex % MathUtility::sqr(BasisFunctionType::nDofsPerBasis())) / BasisFunctionType::nDofsPerBasis());
  dof_no_t localX = dofIndex % BasisFunctionType::nDofsPerBasis();

  return dofsPerPlane * (elementZ * averageNDofsPerElement1D + localZ)
    + dofsPerRow * (elementY * averageNDofsPerElement1D + localY)
    + averageNDofsPerElement1D * elementX + localX;
}

//! get all dofs of a specific node for 3D
template<typename MeshType,typename BasisFunctionType>
void BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const
{
  dofGlobalNos.reserve(dofGlobalNos.size() + BasisOnMeshBaseDim<3,BasisFunctionType>::nDofsPerNode());
  for (int i=0; i<BasisOnMeshBaseDim<3,BasisFunctionType>::nDofsPerNode(); i++)
  {
    dofGlobalNos.push_back(BasisOnMeshBaseDim<3,BasisFunctionType>::nDofsPerNode() * nodeGlobalNo + i);
  }
}

//! get all dofs of a specific node for 3D
template<typename MeshType,typename BasisFunctionType>
void BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getNodeDofs(node_no_t nodeGlobalNo, std::array<dof_no_t,BasisOnMeshBaseDim<3,BasisFunctionType>::nDofsPerNode()> &dofGlobalNos) const
{
  for (int i=0; i<BasisOnMeshBaseDim<3,BasisFunctionType>::nDofsPerNode(); i++)
  {
    dofGlobalNos[i] = BasisOnMeshBaseDim<3,BasisFunctionType>::nDofsPerNode() * nodeGlobalNo + i;
  }
}

//! get the dof no of the specified dof at the node for 3D
template<typename MeshType,typename BasisFunctionType>
dof_no_t BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getNodeDofNo(node_no_t nodeGlobalNo, int dofIndex) const
{
  return BasisOnMeshBaseDim<3,BasisFunctionType>::nDofsPerNode() * nodeGlobalNo + dofIndex;
}

// element-local nodeIndex to local nodeNo for 1D
template<typename MeshType,typename BasisFunctionType>
node_no_t BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getNodeNo(element_no_t elementNo, int nodeIndex) const
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
  return BasisOnMeshFunction<MeshType,BasisFunctionType>::averageNNodesPerElement() * elementNo + nodeIndex;
}

// element-local nodeIndex to local nodeNo for 2D
template<typename MeshType,typename BasisFunctionType>
node_no_t BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
getNodeNo(element_no_t elementNo, int nodeIndex) const
{
  // L linear  quadratic  H cubic
  //           6 7 8
  // 2 3       3 4 5      2 3
  // 0 1       0 1 2      0 1
  // nNodesPerElement:
  // 4         9          4

  // since this implementation is for structured meshes only, the number of elements in each coordinate direction is given

  const std::array<element_no_t, MeshType::dim()> &nElements = this->nElementsPerCoordinateDirection_;
  int averageNNodesPerElement1D = BasisOnMeshBaseDim<1,BasisFunctionType>::averageNNodesPerElement();
  int nNodesPerElement1D = BasisOnMeshBaseDim<1,BasisFunctionType>::nNodesPerElement();
  node_no_t nodesPerRow = (averageNNodesPerElement1D * nElements[0] + 1);
  element_no_t elementX = element_no_t(elementNo % nElements[0]);
  element_no_t elementY = element_no_t(elementNo / nElements[0]);
  dof_no_t localX = nodeIndex % nNodesPerElement1D;
  dof_no_t localY = dof_no_t(nodeIndex / nNodesPerElement1D);

  return nodesPerRow * (elementY * averageNNodesPerElement1D + localY)
    + averageNNodesPerElement1D * elementX + localX;
}

// element-local nodeIndex to local nodeNo for 3D
template<typename MeshType,typename BasisFunctionType>
node_no_t BasisOnMeshNumbers<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getNodeNo(element_no_t elementNo, int nodeIndex) const
{
  // since this implementation is for structured meshes only, the number of elements in each coordinate direction is given

  const std::array<element_no_t, MeshType::dim()> &nElements = this->nElementsPerCoordinateDirection_;
  int averageNNodesPerElement1D = BasisOnMeshBaseDim<1,BasisFunctionType>::averageNNodesPerElement();
  int nNodesPerElement1D = BasisOnMeshBaseDim<1,BasisFunctionType>::nNodesPerElement();
  node_no_t nodesPerRow = (averageNNodesPerElement1D * nElements[0] + 1);
  node_no_t nodesPerPlane = (averageNNodesPerElement1D * nElements[1] + 1) * nodesPerRow;

  element_no_t elementZ = element_no_t(elementNo / (nElements[0] * nElements[1]));
  element_no_t elementY = element_no_t((elementNo % (nElements[0] * nElements[1])) / nElements[0]);
  element_no_t elementX = elementNo % nElements[0];
  dof_no_t localZ = dof_no_t(nodeIndex / MathUtility::sqr(nNodesPerElement1D));
  dof_no_t localY = dof_no_t((nodeIndex % MathUtility::sqr(nNodesPerElement1D)) / nNodesPerElement1D);
  dof_no_t localX = nodeIndex % nNodesPerElement1D;

  return nodesPerPlane * (elementZ * averageNNodesPerElement1D + localZ)
    + nodesPerRow * (elementY * averageNNodesPerElement1D + localY)
    + averageNNodesPerElement1D * elementX + localX;
}

};  // namespace