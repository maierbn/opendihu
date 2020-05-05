#include "function_space/08_function_space_nodes.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"

namespace FunctionSpace
{

template<typename MeshType, typename BasisFunctionType>
std::array<dof_no_t,FunctionSpaceFunction<MeshType,BasisFunctionType>::nNodesPerElement()> FunctionSpaceNodes<MeshType,BasisFunctionType>::
getElementNodeNos(element_no_t elementNo) const
{
  std::array<dof_no_t,FunctionSpaceFunction<MeshType,BasisFunctionType>::nNodesPerElement()> nodes;
  for (int nodeIndex = 0; nodeIndex < FunctionSpaceFunction<MeshType,BasisFunctionType>::nNodesPerElement(); nodeIndex++)
  {
    nodes[nodeIndex] = this->getNodeNo(elementNo, nodeIndex);
  }
  return nodes;
}

//! get the face that is defined by the dofs in the element
template<typename MeshType, typename BasisFunctionType>
Mesh::face_t FunctionSpaceNodes<MeshType,BasisFunctionType>::
getFaceFromElementalDofNos(std::array<int,FunctionSpaceBaseDim<MeshType::dim()-1,BasisFunctionType>::nDofsPerElement()> elementalDofNos)
{
  const int nDofsPerNode = FunctionSpaceBaseDim<MeshType::dim(),BasisFunctionType>::nDofsPerNode();

  int nNodes0Minus = 0;
  int nNodes0Plus = 0;
  int nNodes1Minus = 0;
  int nNodes1Plus = 0;
  int nNodes2Minus = 0;
  int nNodes2Plus = 0;
  
  node_no_t nodeNo = 0;
  int coordinateZ = 0;
  int coordinateY = 0;
  int coordinateX = 0;

  
  switch(MeshType::dim())
  {
  case 3:
    // 3D element
    
    for (int i = 0; i < FunctionSpaceBaseDim<MeshType::dim()-1,BasisFunctionType>::nDofsPerElement(); i++)
    {
      nodeNo = elementalDofNos[i] / nDofsPerNode;
      coordinateZ = nodeNo / 4;
      coordinateY = (nodeNo % 4) / 2;
      coordinateX = nodeNo % 2;

      if (coordinateX == 0)
        nNodes0Minus++;
      if (coordinateX == 1)
        nNodes0Plus++;
      if (coordinateY == 0)
        nNodes1Minus++;
      if (coordinateY == 1)
        nNodes1Plus++;
      if (coordinateZ == 0)
        nNodes2Minus++;
      if (coordinateZ == 1)
        nNodes2Plus++;
    }
    
    if (nNodes0Minus == 4)
      return Mesh::face_t::face0Minus;
    else if (nNodes0Plus == 4)
      return Mesh::face_t::face0Plus;
    else if (nNodes1Minus == 4)
      return Mesh::face_t::face1Minus;
    else if (nNodes1Plus == 4)
      return Mesh::face_t::face1Plus;
    else if (nNodes2Minus == 4)
      return Mesh::face_t::face2Minus;
    else if (nNodes2Plus == 4)
      return Mesh::face_t::face2Plus;

    LOG(FATAL) << "Error in getFaceFromElementalDofNos";

  case 2:

    // 2D element
    
    for (int i = 0; i < FunctionSpaceBaseDim<MeshType::dim()-1,BasisFunctionType>::nDofsPerElement(); i++)
    {
      nodeNo = elementalDofNos[i] / nDofsPerNode;
      coordinateY = nodeNo / 2;
      coordinateX = nodeNo % 2;
      if (coordinateX == 0)
        nNodes0Minus++;
      if (coordinateX == 1)
        nNodes0Plus++;
      if (coordinateY == 0)
        nNodes1Minus++;
      if (coordinateY == 1)
        nNodes1Plus++;
    }

    if (nNodes0Minus == 2)
      return Mesh::face_t::face0Minus;
    else if (nNodes0Plus == 2)
      return Mesh::face_t::face0Plus;
    else if (nNodes1Minus == 2)
      return Mesh::face_t::face1Minus;
    else if (nNodes1Plus == 2)
      return Mesh::face_t::face1Plus;
    
    LOG(FATAL) << "Error in getFaceFromElementalDofNos";

  default:
    // 1D line
    if (elementalDofNos[0] < nDofsPerNode)
    {
      return Mesh::face_t::face0Minus;
    }
    return Mesh::face_t::face0Plus;
  }
}

template<typename MeshType,int D,int order>
dof_no_t FunctionSpaceNodes<MeshType, BasisFunction::CompletePolynomialOfDimensionAndOrder<D,order>>::
nDofsLocal() const
{
  return this->nElementsLocal() * this->nDofsPerElement();
}

template<typename MeshType,int D,int order>
FunctionSpaceNodes<MeshType, BasisFunction::CompletePolynomialOfDimensionAndOrder<D,order>>::
FunctionSpaceNodes(PythonConfig specificSettings, bool noGeometryField) :
  FunctionSpaceFunction<MeshType, BasisFunction::CompletePolynomialOfDimensionAndOrder<MeshType::dim(),order>>::FunctionSpaceFunction(specificSettings)
{
}
} // namespace
