#include "basis_on_mesh/03_basis_on_mesh_dofs.h"

#include <cmath>
#include <array>
#include <string>
#include <cassert>

#include "easylogging++.h"
#include "utility/string_utility.h"
#include "field_variable/field_variable_regular_fixed.h"
#include "field_variable/field_variable_structured_deformable.h"
#include "field_variable/field_variable_unstructured_deformable.h"

namespace BasisOnMesh
{
 
using namespace StringUtility;

// element-local dofIndex to global dofNo for 1D
template<typename MeshType,typename BasisFunctionType>
int BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getDofNo(element_idx_t elementNo, int dofIndex) const
{
  return BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>>::getDofNo(this->nElements_, elementNo, dofIndex);
}

template<typename MeshType,typename BasisFunctionType>
int BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getDofNo(std::array<element_idx_t, MeshType::dim()> nElements, element_idx_t elementNo, int dofIndex)
{
  // L linear  L quadratic  H cubic
  // 0 1       0 1 2        0,1 2,3
  // averageNDofsPerElement: 
  // 1         2            2
  return BasisOnMeshFunction<MeshType,BasisFunctionType>::averageNDofsPerElement() * elementNo + dofIndex;
}

//! get all dofs of a specific node for 1D
template<typename MeshType,typename BasisFunctionType>
void BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getNodeDofs(node_idx_t nodeGlobalNo, std::vector<int> dofGlobalNos) const
{
  for (int i=0; i<BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerNode(); i++)
  {
    dofGlobalNos.push_back(BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerNode() * nodeGlobalNo + i);
  }
}

// element-local dofIndex to global dofNo for 2D
template<typename MeshType,typename BasisFunctionType>
int BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
getDofNo(element_idx_t elementNo, int dofIndex) const
{
  return BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>>::getDofNo(this->nElements_, elementNo, dofIndex);
}

// element-local dofIndex to global dofNo for 2D
template<typename MeshType,typename BasisFunctionType>
int BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
getDofNo(std::array<element_idx_t, MeshType::dim()> nElements, element_idx_t elementNo, int dofIndex)
{
  // L linear  quadratic  H cubic
  //           6 7 8
  // 2 3       3 4 5      45 67
  // 0 1       0 1 2      01 23
  // nDofsPerBasis:
  // 2         3          4
  // averageNDofsPerElement:
  // 1         4          2
  
  int averageNDofsPerElement1D = BasisOnMeshBaseDim<1,BasisFunctionType>::averageNDofsPerElement();
  int dofsPerRow = (averageNDofsPerElement1D * nElements[0] + BasisFunctionType::nDofsPerNode());
  int elementX = int(elementNo % nElements[0]);
  int elementY = int(elementNo / nElements[0]);
  int localX = dofIndex % BasisFunctionType::nDofsPerBasis();
  int localY = int(dofIndex / BasisFunctionType::nDofsPerBasis());
  
  VLOG(2) << "  dof " << elementNo << ":" << dofIndex << ", element: ("<<elementX<<","<<elementY<<"), dofsPerRow="<<dofsPerRow<<", local: ("<<localX<<","<<localY<<")";
  
  return dofsPerRow * (elementY * averageNDofsPerElement1D + localY) 
    + averageNDofsPerElement1D * elementX + localX;
}

//! get all dofs of a specific node for 2D
template<typename MeshType,typename BasisFunctionType>
void BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
getNodeDofs(node_idx_t nodeGlobalNo, std::vector<int> dofGlobalNos) const
{
  for (int i=0; i<BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerNode(); i++)
  {
    dofGlobalNos.push_back(BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerNode() * nodeGlobalNo + i);
  }
}
  
// element-local dofIndex to global dofNo for 3D
template<typename MeshType,typename BasisFunctionType>
int BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getDofNo(element_idx_t elementNo, int dofIndex) const
{
  return BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>>::getDofNo(this->nElements_, elementNo, dofIndex);
}
 
// element-local dofIndex to global dofNo for 3D
template<typename MeshType,typename BasisFunctionType>
int BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getDofNo(std::array<element_idx_t, MeshType::dim()> nElements, element_idx_t elementNo, int dofIndex)
{
  // L linear  quadratic  H cubic
  //           6 7 8
  // 2 3       3 4 5      45 67
  // 0 1       0 1 2      01 23
  // nDofsPerBasis:
  // 2         3          4
  // averageNDofsPerElement:
  // 1         4          2
 
  int averageNDofsPerElement1D = BasisOnMeshBaseDim<1,BasisFunctionType>::averageNDofsPerElement();
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

//! get all dofs of a specific node for 3D
template<typename MeshType,typename BasisFunctionType>
void BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getNodeDofs(node_idx_t nodeGlobalNo, std::vector<int> dofGlobalNos) const
{
  for (int i=0; i<BasisOnMeshBaseDim<3,BasisFunctionType>::nDofsPerNode(); i++)
  {
    dofGlobalNos.push_back(BasisOnMeshBaseDim<3,BasisFunctionType>::nDofsPerNode() * nodeGlobalNo + i);
  }
}

// element-local nodeIndex to global nodeNo for 1D
template<typename MeshType,typename BasisFunctionType>
int BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<1,MeshType>> ::
getNodeNo(element_idx_t elementNo, int nodeIndex) const
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

// element-local nodeIndex to global nodeNo for 2D
template<typename MeshType,typename BasisFunctionType>
int BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<2,MeshType>> ::
getNodeNo(element_idx_t elementNo, int nodeIndex) const
{
  // L linear  quadratic  H cubic
  //           6 7 8
  // 2 3       3 4 5      2 3
  // 0 1       0 1 2      0 1
  // nNodesPerElement:
  // 4         9          4
  
  // since this implementation is for structured meshes only, the number of elements in each coordinate direction is given
  const std::array<int,2> nElements{this->nElements(0), this->nElements(1)};
  
  int averageNNodesPerElement1D = BasisOnMeshBaseDim<1,BasisFunctionType>::averageNNodesPerElement();
  int nNodesPerElement1D = BasisOnMeshBaseDim<1,BasisFunctionType>::nNodesPerElement();
  int nodesPerRow = (averageNNodesPerElement1D * nElements[0] + 1);
  int elementX = int(elementNo % nElements[0]);
  int elementY = int(elementNo / nElements[0]);
  int localX = nodeIndex % nNodesPerElement1D;
  int localY = int(nodeIndex / nNodesPerElement1D);
  
  return nodesPerRow * (elementY * averageNNodesPerElement1D + localY) 
    + averageNNodesPerElement1D * elementX + localX;
}

// element-local nodeIndex to global nodeNo for 3D
template<typename MeshType,typename BasisFunctionType>
int BasisOnMeshDofs<MeshType,BasisFunctionType,Mesh::isStructuredWithDim<3,MeshType>> ::
getNodeNo(element_idx_t elementNo, int nodeIndex) const
{
  // since this implementation is for structured meshes only, the number of elements in each coordinate direction is given
  const std::array<int,3> nElements{this->nElements(0), this->nElements(1), this->nElements(2)};
  
  int averageNNodesPerElement1D = BasisOnMeshBaseDim<1,BasisFunctionType>::averageNNodesPerElement();
  int nNodesPerElement1D = BasisOnMeshBaseDim<1,BasisFunctionType>::nNodesPerElement();
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


template<int D,typename BasisFunctionType>
element_idx_t BasisOnMeshDofs<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::
nElements() const
{
  return this->nElements_;
}

template<int D,typename BasisFunctionType>
BasisOnMeshDofs<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::
BasisOnMeshDofs(PyObject *settings) :
  BasisOnMeshJacobian<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::BasisOnMeshJacobian(settings)
{
  std::string filenameExelem = PythonUtility::getOptionString(settings, "exelem", "input.exelem");
  std::string filenameExnode = PythonUtility::getOptionString(settings, "exnode", "input.exnode");
  
  LOG(TRACE) << "BasisOnMeshDofs constructor";
  
  // ------------------------------------------------------------------------
  // read in exelem file
  this->parseExelemFile(filenameExelem);
  
  // ------------------------------------------------------------------------
  // read in exnode file
  this->parseExnodeFile(filenameExnode);
  
  // remap names of field variables if specified in config
  this->remapFieldVariables(settings);
  
  
  // eliminate scale factors
  //this->eliminateScaleFactors();
}


template<int D,typename BasisFunctionType>
void BasisOnMeshDofs<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::
initialize()
{
  LOG(DEBUG) << "   retrieve this pointer ";
  std::shared_ptr<BasisOnMeshDofs<Mesh::UnstructuredDeformable<D>,BasisFunctionType>> ptr = this->shared_from_this();
  
  assert(ptr != nullptr);
  
  LOG(DEBUG) << "   cast this pointer ";
  std::shared_ptr<BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>> self = std::static_pointer_cast<BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>(ptr);
  
  assert(self != nullptr);
  assert(fieldVariable_.find("geometry") != fieldVariable_.end());
  
  this->fieldVariable_.at("geometry")->setMesh(self);
}

template<int D,typename BasisFunctionType>
void BasisOnMeshDofs<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::
remapFieldVariables(PyObject *settings)
{
  // remap name of field variables to different names if specified
  if (PythonUtility::containsKey(settings, "remap"))
  {
    std::string keyString = "remap";
    std::pair<std::string, std::string> dictItem 
      = PythonUtility::getOptionDictBegin<std::string, std::string>(settings, keyString);
    
    for (; !PythonUtility::getOptionDictEnd(settings, keyString); 
        PythonUtility::getOptionDictNext<std::string, std::string>(settings, keyString, dictItem))
    {
      std::string key = dictItem.first;
      std::string value = dictItem.second;
          
      if (this->fieldVariable_.find(key) != this->fieldVariable_.end())
      {
        std::shared_ptr<FieldVariableType> fieldVariable = this->fieldVariable_[key];
        this->fieldVariable_.erase(key);
        this->fieldVariable_[value] = fieldVariable;
      }
      else
      {
        LOG(WARNING) << "Remap of field variable from \"" << key << "\" to \"" << value << "\" failed: no such entry.";
      }
    }
  }
  
  // if there is no field with name "geometry"
  if (this->fieldVariable_.find("geometry") == this->fieldVariable_.end())
  { 
    bool geometryFieldFound = false;
    // search for a geometry field 
    for(auto &fieldVariableEntry : this->fieldVariable_)
    {
      if (fieldVariableEntry.second->isGeometryField())
      {
        LOG(WARNING) << "Remap geometry field variable from \"" << fieldVariableEntry.first << "\" to \"geometry\".";
      
        std::shared_ptr<FieldVariableType> fieldVariable = fieldVariableEntry.second;
        this->fieldVariable_.erase(fieldVariableEntry.first);
        this->fieldVariable_["geometry"] = fieldVariable;
        
        geometryFieldFound = true;
        break;
      }
    }
    
    // output field variables
    for(auto &fieldVariableEntry : this->fieldVariable_)
    {
      VLOG(1) << *fieldVariableEntry.second;
    }
    
    if (!geometryFieldFound)
    {
      LOG(FATAL) << "The specified Exfiles contain no geometry field. The field must be named \"geometry\" or have the type \"coordinates\".";
    }
  }
}

template<int D,typename BasisFunctionType>
int BasisOnMeshDofs<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::
getDofNo(element_idx_t elementNo, int dofIndex) const
{
  if (this->fieldVariable_.find("geometry") == this->fieldVariable_.end())
    LOG(FATAL) << "Mesh contains no field variable \"geometry\". Use remap to create one!";
  
  return this->fieldVariable_.at("geometry")->getDofNo(elementNo, dofIndex);
}

template<int D,typename BasisFunctionType>
int BasisOnMeshDofs<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::
getNodeNo(element_idx_t elementNo, int nodeIndex) const
{
  return this->elementToNodeMapping_->getElement(elementNo).nodeGlobalNo[nodeIndex];
}

//! get all dofs of a specific node, unstructured mesh
template<int D,typename BasisFunctionType>
void BasisOnMeshDofs<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::
getNodeDofs(node_idx_t nodeGlobalNo, std::vector<int> dofGlobalNos) const
{
  if (this->fieldVariable_.find("geometry") == this->fieldVariable_.end())
    LOG(FATAL) << "Mesh contains no field variable \"geometry\". Use remap to create one!";
  
  std::vector<int> &nodeDofs = dofGlobalNos.insert(dofGlobalNos.end(), 
    this->fieldVariable_.at("geometry")->nodeToDofMapping()->getNodeDofs(nodeGlobalNo));
  dofGlobalNos.reserve(dofGlobalNos.size() + nodeDofs.size());
  
  for(int dof : nodeDofs)
  {
    dofGlobalNos.push_back(dof);
  }
}

template<int D,typename BasisFunctionType>
void BasisOnMeshDofs<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::
eliminateScaleFactors()
{
  // loop over field variables
  for (auto fieldVariable : this->fieldVariable_)
  {
    fieldVariable.second->eliminateScaleFactors();
  }
}

template<int D,typename BasisFunctionType>
int BasisOnMeshDofs<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::
nDofs() const
{
  return nDofs_;
}

template<int D,typename BasisFunctionType>
void BasisOnMeshDofs<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::
addNonGeometryFieldVariables(std::vector<std::shared_ptr<FieldVariableType>> &fieldVariables)
{
  // loop over field variables
  for (auto fieldVariable : this->fieldVariable_)
  {
    if (!fieldVariable.second->isGeometryField())
      fieldVariables.push_back(fieldVariable.second);
  }
}
  
};  // namespace
