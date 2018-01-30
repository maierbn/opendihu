#include "field_variable/field_variable_regular_fixed.h"

#include <sstream>
#include "easylogging++.h"
#include "utility/string_utility.h"

namespace FieldVariable
{
 
 
template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::
setMeshWidth(std::array<double, D> &meshWidth)
{
  meshWidth_ = meshWidth;
}


template<int D, typename BasisFunctionType>
double FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::
meshWidth(int dimension) const
{
  return this->meshWidth_[dimension];
}

//! for a specific component, get values from their global dof no.s
template<int D,typename BasisFunctionType>
template<int N>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::
getValues(std::string component, std::array<int,N> dofGlobalNo, std::array<double,N> &values)
{
  if (!this->isGeometryField_)
  {
    FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::getValues(component, dofGlobalNo, values);
    return;
  }
  
  // for geometry field compute information
  const int componentIndex = this->componentIndex_[component];
  
  const node_idx_t nNodesInXDirection = this->mesh_->nNodes(0);
  const node_idx_t nNodesInYDirection = this->mesh_->nNodes(1);
  const int nDofsPerNode = BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>::nDofsPerNode();
  
 
  // loop over entries in values to be filled
  for (int i=0; i<N; i++)
  {
    int nodeNo = int(dofGlobalNo[i] / nDofsPerNode);
    int nodeLocalDofIndex = int(dofGlobalNo[i] % nDofsPerNode);
    
    if (nodeLocalDofIndex > 0)   // if this is a derivative of Hermite, set to 0
    {
      values[i] = 0;
    }
    else
    {
      if (componentIndex == 0)   // x direction
      {
        values[i] = (nodeNo % nNodesInXDirection) * this->meshWidth_[0];
      }
      else if (componentIndex == 1)   // y direction
      {
        values[i] = (int(nodeNo / nNodesInXDirection) % nNodesInYDirection) * this->meshWidth_[1];
      }
      else if (componentIndex == 2)   // z direction
      {
        values[i] = int(nodeNo / (nNodesInXDirection*nNodesInYDirection)) * this->meshWidth_[2];
      }
    }
  }
}

/*
//! get values from their global dof no.s for all components
template<int D,typename BasisFunctionType>
template<int N, int nComponents>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::
getValues(std::array<int,N> dofGlobalNo, std::array<std::array<double,nComponents>,N> &values)
{
  if (!this->isGeometryField_)
  {
    FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::getValues(dofGlobalNo, values);
    return;
  }

  const node_idx_t nNodesInXDirection = this->mesh_->nNodes(0);
  const node_idx_t nNodesInYDirection = this->mesh_->nNodes(1);
  const int nDofsPerNode = BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>::nDofsPerNode();
  
  // loop over entries in values to be filled
  for (int i=0; i<N; i++)
  {
    int nodeNo = int(dofGlobalNo[i] / nDofsPerNode);
    int nodeLocalDofIndex = int(dofGlobalNo[i] % nDofsPerNode);
    
    if (nodeLocalDofIndex > 0)   // if this is a derivative of Hermite, set to 0
    {
      values[i][0] = 0;
      values[i][1] = 0;
      values[i][2] = 0;
    }
    else 
    {
      // x direction
      values[i][0] = (nodeNo % nNodesInXDirection) * this->meshWidth_[0];
      
      // y direction
      values[i][1] = (int(nodeNo / nNodesInXDirection) % nNodesInYDirection) * this->meshWidth_[1];
      
      // z direction
      values[i][2] = int(nodeNo / (nNodesInXDirection*nNodesInYDirection)) * this->meshWidth_[2];
    }
  }
}
*/

//! for a specific component, get the values corresponding to all element-local dofs
template<int D,typename BasisFunctionType>
template<int N>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::
getElementValues(std::string component, element_idx_t elementNo, 
                 std::array<double,BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>::nDofsPerElement()> &values)
{
  if (!this->isGeometryField_)
  {
    FieldVariableStructured<BasisOnMeshType>::getValues(component, elementNo, values);
    return;
  }
  
  int componentIndex = this->componentIndex_[component];
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  
  // get the global dof no.s of the element
  std::array<int,nDofsPerElement> dofGlobalNo;
  for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
  {
    dofGlobalNo[dofIndex] = BasisOnMeshType::dofNo(elementNo, dofIndex);
  }
  
  this->getValues<nDofsPerElement,D,BasisFunctionType>(component, dofGlobalNo, values);
}
/*
//! get the values corresponding to all element-local dofs for all components
template<int D,typename BasisFunctionType>
template<int nComponents>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::
getElementValues(element_idx_t elementNo, std::array<std::array<double,nComponents>,BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>::nDofsPerElement()> &values)
{
  if (!this->isGeometryField_)
  {
    FieldVariableStructured<BasisOnMeshType>::getValues(elementNo, values);
    return;
  }
  
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  
  // get the global dof no.s of the element
  std::array<int,nDofsPerElement> dofGlobalNo;
  for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
  {
    dofGlobalNo[dofIndex] = BasisOnMeshType::dofNo(elementNo, dofIndex);
  }
  
  this->getValues<nDofsPerElement,nComponents,D,BasisFunctionType>(dofGlobalNo, values);
}
*/
//! for a specific component, get a single value from global dof no.
template<int D,typename BasisFunctionType>
double FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::
getValue(std::string component, node_idx_t dofGlobalNo)
{
  if (!this->isGeometryField_)
  {
    return FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::getValue(component, dofGlobalNo);
  }
  
  // for geometry field compute information
  const int componentIndex = this->componentIndex_[component];
  
  const node_idx_t nNodesInXDirection = this->mesh_->nNodes(0);
  const node_idx_t nNodesInYDirection = this->mesh_->nNodes(1);
  const int nDofsPerNode = BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>::nDofsPerNode();
  
  double value = 0;
  int nodeNo = int(dofGlobalNo / nDofsPerNode);
  int nodeLocalDofIndex = int(dofGlobalNo % nDofsPerNode);
  
  if (nodeLocalDofIndex == 0)   // if this is not a derivative of Hermite (in which case it would be set to 0)
  {
    if (componentIndex == 0)   // x direction
    {
      value = (nodeNo % nNodesInXDirection) * this->meshWidth_[0];
    }
    else if (componentIndex == 1)   // y direction
    {
      value = (int(nodeNo / nNodesInXDirection) % nNodesInYDirection) * this->meshWidth_[1];
    }
    else if (componentIndex == 2)   // z direction
    {
      value = int(nodeNo / (nNodesInXDirection*nNodesInYDirection)) * this->meshWidth_[2];
    }
  }
  return value;
}

//! get a single value from global dof no. for all components
template<int D,typename BasisFunctionType>
template<int nComponents>
std::array<double,nComponents> FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::
getValue(node_idx_t dofGlobalNo)
{
  if (!this->isGeometryField_)
  {
    return FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::template getValue<nComponents>(dofGlobalNo);
  }
  
  const node_idx_t nNodesInXDirection = this->mesh_->nNodes(0);
  const node_idx_t nNodesInYDirection = this->mesh_->nNodes(1);
  const int nDofsPerNode = BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>::nDofsPerNode();
  
  int nodeNo = int(dofGlobalNo / nDofsPerNode);
  int nodeLocalDofIndex = int(dofGlobalNo % nDofsPerNode);
  
  std::array<double,nComponents> value;
  if (nodeLocalDofIndex > 0)   // if this is a derivative of Hermite, set to 0
  {
    value[0] = 0;
    value[1] = 0;
    value[2] = 0;
  }
  else 
  {
    // x direction
    value[0] = (nodeNo % nNodesInXDirection) * this->meshWidth_[0];
    
    // y direction
    value[1] = (int(nodeNo / nNodesInXDirection) % nNodesInYDirection) * this->meshWidth_[1];
    
    // z direction
    value[2] = int(nodeNo / (nNodesInXDirection*nNodesInYDirection)) * this->meshWidth_[2];
  }
  return value;
}

//! write a exelem file header to a stream, for a particular element
template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::
outputHeaderExelem(std::ostream &file, element_idx_t currentElementGlobalNo, int fieldVariableNo)
{
  // use the implementation of FieldVariableStructured
  FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::
    outputHeaderExelem(file, currentElementGlobalNo, fieldVariableNo);
}

//! write a exelem file header to a stream, for a particular element
template<int D, typename BasisFunctionType>
void FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::
outputHeaderExnode(std::ostream &file, node_idx_t currentNodeGlobalNo, int &valueIndex, int fieldVariableNo)
{
  // use the implementation of FieldVariableStructured
  FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::
    outputHeaderExnode(file, currentNodeGlobalNo, valueIndex, fieldVariableNo);
}

//! tell if 2 elements have the same exfile representation, i.e. same number of versions
template<int D, typename BasisFunctionType>
bool FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::
haveSameExfileRepresentation(element_idx_t element1, element_idx_t element2)
{
  // use the implementation of FieldVariableStructured
  return FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::
    haveSameExfileRepresentation(element1, element2);
}

//! get the internal PETSc vector values. The meaning of the values is instance-dependent (different for different BasisOnMeshTypes)
template<int D, typename BasisFunctionType>
Vec &FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::
values()
{
  // use the implementation of FieldVariableStructured
  return FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::values();
}

//! get the number of components
template<int D, typename BasisFunctionType>
int FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::
nComponents() const
{
  // use the implementation of FieldVariableStructured
  return FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::nComponents();
}

//! get the number of elements
template<int D, typename BasisFunctionType>
std::array<int, BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>::Mesh::dim()> FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::
nElementsPerDimension() const
{
  // use the implementation of FieldVariableStructured
  return FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::nElementsPerDimension();
}

//! get the number of elements
template<int D, typename BasisFunctionType>
element_idx_t FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::
nElements() const
{
  return this->nElements();
}

//! get the names of the components that are part of this field variable
template<int D, typename BasisFunctionType>
std::vector<std::string> FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::
componentNames() const
{
  // use the implementation of FieldVariableStructured
  return FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::componentNames();
}
};
