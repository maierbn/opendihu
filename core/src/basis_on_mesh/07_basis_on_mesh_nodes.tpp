#include "basis_on_mesh/07_basis_on_mesh_nodes.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"

namespace BasisOnMesh
{

template<typename MeshType, typename BasisFunctionType>
std::array<dof_no_t,BasisOnMeshFunction<MeshType,BasisFunctionType>::nNodesPerElement()> BasisOnMeshNodes<MeshType,BasisFunctionType>::
getElementNodeNos(element_no_t elementNo) const
{
  std::array<dof_no_t,BasisOnMeshFunction<MeshType,BasisFunctionType>::nNodesPerElement()> nodes;
  for (int nodeIndex = 0; nodeIndex < BasisOnMeshFunction<MeshType,BasisFunctionType>::nNodesPerElement(); nodeIndex++)
  {
    nodes[nodeIndex] = this->getNodeNo(elementNo, nodeIndex);
  }
  return nodes;
}

template<typename MeshType,int D,int order>
dof_no_t BasisOnMeshNodes<MeshType, BasisFunction::CompletePolynomialOfDimensionAndOrder<D,order>>::
nDofs() const
{
  return this->nElements() * this->nDofsPerElement();
}

template<typename MeshType,int D,int order>
BasisOnMeshNodes<MeshType, BasisFunction::CompletePolynomialOfDimensionAndOrder<D,order>>::
BasisOnMeshNodes(PyObject *specificSettings, bool noGeometryField) :
  BasisOnMeshFunction<MeshType, BasisFunction::CompletePolynomialOfDimensionAndOrder<MeshType::dim(),order>>::BasisOnMeshFunction(specificSettings)
{
}
};  // namespace