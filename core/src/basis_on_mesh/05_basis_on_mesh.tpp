#include "basis_on_mesh/05_basis_on_mesh.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"

namespace BasisOnMesh
{

template<typename MeshType, typename BasisFunctionType>
std::array<int,BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement()> BasisOnMesh<MeshType,BasisFunctionType>::
getElementDofNos(element_idx_t elementNo) const
{
  std::array<int,BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement()> dof;
  for (int dofIndex = 0; dofIndex < BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement(); dofIndex++)
  {
    dof[dofIndex] = this->getDofNo(elementNo, dofIndex);
  }
  return dof;
}

template<typename MeshType, typename BasisFunctionType>
std::array<int,BasisOnMeshFunction<MeshType,BasisFunctionType>::nNodesPerElement()> BasisOnMesh<MeshType,BasisFunctionType>::
getElementNodeNos(element_idx_t elementNo) const
{
  std::array<int,BasisOnMeshFunction<MeshType,BasisFunctionType>::nNodesPerElement()> nodes;
  for (int nodeIndex = 0; nodeIndex < BasisOnMeshFunction<MeshType,BasisFunctionType>::nNodesPerElement(); nodeIndex++)
  {
    nodes[nodeIndex] = this->getNodeNo(elementNo, nodeIndex);
  }
  return nodes;
}

template<typename MeshType, typename BasisFunctionType>
std::array<std::array<double,MeshType::dim()>,BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement()> BasisOnMesh<MeshType,BasisFunctionType>::
getGradPhi(std::array<double,MeshType::dim()> xi) const
{
  std::array<std::array<double,MeshType::dim()>,BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement()> gradPhi;
  for (int dofIndex = 0; dofIndex < BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement(); dofIndex++)
  {
    gradPhi[dofIndex] = this->gradPhi(dofIndex, xi);
  }
  return gradPhi;
}

//! create a non-geometry field field variable with no values being set, with given component names
template<typename MeshType, typename BasisFunctionType>
std::shared_ptr<FieldVariable::FieldVariable<BasisOnMesh<MeshType,BasisFunctionType>>> BasisOnMesh<MeshType,BasisFunctionType>::
createFieldVariable(std::string name, std::vector<std::string> componentNames)
{
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMesh<MeshType,BasisFunctionType>>> fieldVariable
    = std::make_shared<FieldVariable::FieldVariable<BasisOnMesh<MeshType,BasisFunctionType>>>();
  fieldVariable->initializeFromFieldVariable(this->geometryField(), name, componentNames);
  
  return fieldVariable;
}

//! create a non-geometry field field variable with no values being set, with given number of components, the component names will be the numbers
template<typename MeshType, typename BasisFunctionType>
std::shared_ptr<FieldVariable::FieldVariable<BasisOnMesh<MeshType,BasisFunctionType>>> BasisOnMesh<MeshType,BasisFunctionType>::
createFieldVariable(std::string name, int nComponents)
{
  std::vector<std::string> componentNames(nComponents);
  for (int i=0; i<nComponents; i++)
  {
    componentNames[i] = std::to_string(i);
  }
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMesh<MeshType,BasisFunctionType>>> fieldVariable 
    = this->createFieldVariable(name, componentNames); 
  
  return fieldVariable;
}

};  // namespace