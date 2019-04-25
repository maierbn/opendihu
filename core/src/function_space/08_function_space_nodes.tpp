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
