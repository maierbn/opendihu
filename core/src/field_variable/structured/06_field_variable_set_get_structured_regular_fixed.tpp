#include "field_variable/structured/06_field_variable_set_get_structured_regular_fixed.h"

#include <sstream>
#include "easylogging++.h"
#include "utility/string_utility.h"

namespace FieldVariable
{

//! for a specific component, get all values
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetRegularFixed<BasisOnMeshType,nComponents>::
getValues(int componentNo, std::vector<double> &values, bool onlyNodalValues)
{
  // if this is not a geometry field get the stored values
  if (!this->isGeometryField_)
  {
    FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
      getValues(componentNo, values, onlyNodalValues);
    return;
  }
  
  // for geometry field compute information
  node_no_t nNodesInXDirection = this->mesh_->nNodes(0);
  node_no_t nNodesInYDirection = this->mesh_->nNodes(1);
  node_no_t nNodesInZDirection = this->mesh_->nNodes(2);
  
  const int D = BasisOnMeshType::dim();
  if (D < 2)
    nNodesInYDirection = 1;
  if (D < 3)
    nNodesInZDirection = 1;
  
  const int nDofsPerNode = BasisOnMeshType::nDofsPerNode();
 
  // determine the number of values to be retrived which is half the number of dofs for Hermite with only nodal values
  dof_no_t nValues = this->mesh_->nDofs();
  if (onlyNodalValues)
    // if the basis function is Hermite
    if (std::is_same<typename BasisOnMeshType::BasisFunction, BasisFunction::Hermite>::value)
      nValues = this->mesh_->nDofs() / 2;
 
  LOG(DEBUG) << "getValues, n dofs: " << this->mesh_->nDofs() << ", nValues: " << nValues
    << ", nNodes: " << nNodesInXDirection<<","<<nNodesInYDirection<<","<<nNodesInZDirection
    << ", nDofsPerNode: " << nDofsPerNode;
  
  values.resize(nValues);
  std::size_t vectorIndex = 0;
  
  // loop over all nodes
  for (int nodeZ = 0; nodeZ < nNodesInZDirection; nodeZ++)
  {
    for (int nodeY = 0; nodeY < nNodesInYDirection; nodeY++)
    {
      for (int nodeX = 0; nodeX < nNodesInXDirection; nodeX++)
      {
        if (componentNo == 0)  // "x"
        {
          assert(vectorIndex < values.size());
          values[vectorIndex++] = nodeX * this->meshWidth_;
        }
        else if (componentNo == 1)  // "y"
        {
          assert(vectorIndex < values.size());
          values[vectorIndex++] = nodeY * this->meshWidth_;
        }
        else  // "z"
        {
          assert(vectorIndex < values.size());
          values[vectorIndex++] = nodeZ * this->meshWidth_;
        }
       
        // set derivative of Hermite to 0 for geometry field
        if (!onlyNodalValues)
        {
          for (int dofIndex = 1; dofIndex < nDofsPerNode; dofIndex++)
          {
            values[vectorIndex++] = 0;
          }
        }
      }
    }
  }
}

//! for a specific component, get values from their global dof no.s
template<typename BasisOnMeshType, int nComponents>
template<int N>
void FieldVariableSetGetRegularFixed<BasisOnMeshType,nComponents>::
getValues(int componentNo, std::array<dof_no_t,N> dofGlobalNo, std::array<double,N> &values)
{
  // if this is not a geometry field get the stored values
  if (!this->isGeometryField_)
  {
    FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
      template getValues<N>(componentNo, dofGlobalNo, values);
    return;
  }
  
  // for geometry field compute information
  const node_no_t nNodesInXDirection = this->mesh_->nNodes(0);
  const node_no_t nNodesInYDirection = this->mesh_->nNodes(1);
  const int nDofsPerNode = BasisOnMeshType::nDofsPerNode();
  
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
      if (componentNo == 0)   // x direction
      {
        values[i] = (nodeNo % nNodesInXDirection) * this->meshWidth_;
      }
      else if (componentNo == 1)   // y direction
      {
        values[i] = (int(nodeNo / nNodesInXDirection) % nNodesInYDirection) * this->meshWidth_;
      }
      else  // z direction
      {
        values[i] = int(nodeNo / (nNodesInXDirection*nNodesInYDirection)) * this->meshWidth_;
      }
    }
  }
}

//! get values from their global dof no.s for all components
template<typename BasisOnMeshType, int nComponents>
template<int N>
void FieldVariableSetGetRegularFixed<BasisOnMeshType,nComponents>::
getValues(std::array<dof_no_t,N> dofGlobalNo, std::array<std::array<double,nComponents>,N> &values)
{
  // if this is not a geometry field get the stored values
  if (!this->isGeometryField_)
  {
    FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
      template getValues<N>(dofGlobalNo, values);
    return;
  }

  // for geometry field compute the entries
  const node_no_t nNodesInXDirection = this->mesh_->nNodes(0);
  const node_no_t nNodesInYDirection = this->mesh_->nNodes(1);
  const int nDofsPerNode = BasisOnMeshType::nDofsPerNode();
  
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
      values[i][0] = (nodeNo % nNodesInXDirection) * this->meshWidth_;
      
      // y direction
      values[i][1] = (int(nodeNo / nNodesInXDirection) % nNodesInYDirection) * this->meshWidth_;
      
      // z direction
      values[i][2] = int(nodeNo / (nNodesInXDirection*nNodesInYDirection)) * this->meshWidth_;
    }
  }
}
  
//! for a specific component, get the values corresponding to all element-local dofs 
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetRegularFixed<BasisOnMeshType,nComponents>::
getElementValues(int componentNo, element_no_t elementNo, std::array<double,BasisOnMeshType::nDofsPerElement()> &values)
{
  // if this is not a geometry field get the stored values
  if (!this->isGeometryField_)
  {
    FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
      getElementValues(componentNo, elementNo, values);
    return;
  }
  
  // if this is a geometry field, compute the entries
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  
  // get the element-local dofs of the element
  std::array<dof_no_t,nDofsPerElement> elementDofs = this->mesh_->getElementDofNos(elementNo);
  
  // get the values
  this->getValues<nDofsPerElement>(componentNo, elementDofs, values);
}

//! get the values corresponding to all element-local dofs for all components
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetRegularFixed<BasisOnMeshType,nComponents>::
getElementValues(element_no_t elementNo, std::array<std::array<double,nComponents>,BasisOnMeshType::nDofsPerElement()> &values)
{
  // if this is not a geometry field get the stored values
  if (!this->isGeometryField_)
  {
    FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
      getElementValues(elementNo, values);
    return;
  }
  
  // for geometry field compute the requested values
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
   
  // get the element-local dofs of the element
  std::array<dof_no_t,nDofsPerElement> elementDofs = this->mesh_->getElementDofNos(elementNo);
  
  // compute the corresponding geometry values
  this->getValues<nDofsPerElement>(elementDofs, values);
}
 
//! for a specific component, get a single value from global dof no.
template<typename BasisOnMeshType, int nComponents>
double FieldVariableSetGetRegularFixed<BasisOnMeshType,nComponents>::
getValue(int componentNo, node_no_t dofGlobalNo)
{
  if (!this->isGeometryField_)
  {
    return FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
      getValue(componentNo, dofGlobalNo);
  }
  
  // for geometry field compute information
  const node_no_t nNodesInXDirection = this->mesh_->nNodes(0);
  const node_no_t nNodesInYDirection = this->mesh_->nNodes(1);
  const int nDofsPerNode = BasisOnMeshType::nDofsPerNode();
  
  double value = 0;
  int nodeNo = int(dofGlobalNo / nDofsPerNode);
  int nodeLocalDofIndex = int(dofGlobalNo % nDofsPerNode);
  
  if (nodeLocalDofIndex == 0)   // if this is not a derivative of Hermite (in which case it would be set to 0)
  {
    if (componentNo == 0)   // x direction
    {
      value = (nodeNo % nNodesInXDirection) * this->meshWidth_;
    }
    else if (componentNo == 1)   // y direction
    {
      value = (int(nodeNo / nNodesInXDirection) % nNodesInYDirection) * this->meshWidth_;
    }
    else     // z direction
    {
      value = int(nodeNo / (nNodesInXDirection*nNodesInYDirection)) * this->meshWidth_;
    }
  }
  return value;
}

//! copy the values from another field variable of the same type
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetRegularFixed<BasisOnMeshType,nComponents>::
setValues(FieldVariable<BasisOnMeshType,nComponents> &rhs)
{
  VecCopy(rhs.values_, this->values_);
}

};
