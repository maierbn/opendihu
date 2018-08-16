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
    VLOG(3) << "getValues(componentNo=" << componentNo << "): field variable is not a geometry field, retrieve stored values";
    FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
      getValues(componentNo, values, onlyNodalValues);
    return;
  }

  // for geometry field compute information
  node_no_t nLocalNodesInXDirection = this->mesh_->nNodesLocalWithGhosts(0);
  node_no_t nLocalNodesInYDirection = this->mesh_->nNodesLocalWithGhosts(1);
  node_no_t nLocalNodesInZDirection = this->mesh_->nNodesLocalWithGhosts(2);

  const int D = BasisOnMeshType::dim();
  if (D < 2)
    nLocalNodesInYDirection = 1;
  if (D < 3)
    nLocalNodesInZDirection = 1;

  const int nDofsPerNode = BasisOnMeshType::nDofsPerNode();

  // determine the number of values to be retrived which is lower than the number of dofs for Hermite with only nodal values
  dof_no_t nValues = this->mesh_->nLocalDofs();
  if (onlyNodalValues)
    // if the basis function is Hermite
    if (std::is_same<typename BasisOnMeshType::BasisFunction, BasisFunction::Hermite>::value)
      nValues = this->mesh_->nLocalDofs() / BasisOnMeshType::nDofsPerNode();

  LOG(DEBUG) << "getValues, n dofs: " << this->mesh_->nLocalDofs() << ", nValues: " << nValues
    << ", nNodes: " << nLocalNodesInXDirection<<","<<nLocalNodesInYDirection<<","<<nLocalNodesInZDirection
    << ", nDofsPerNode: " << nDofsPerNode;

  values.resize(nValues);
  std::size_t vectorIndex = 0;

  // loop over all nodes
  for (int nodeZ = 0; nodeZ < nLocalNodesInZDirection; nodeZ++)
  {
    for (int nodeY = 0; nodeY < nLocalNodesInYDirection; nodeY++)
    {
      for (int nodeX = 0; nodeX < nLocalNodesInXDirection; nodeX++)
      {
        if (componentNo == 0)  // "x"
        {
          assert(vectorIndex < values.size());
          values[vectorIndex++] = this->mesh_->meshPartition()->beginNodeGlobal(0) * this->meshWidth_ 
            + nodeX * this->meshWidth_;
        }
        else if (componentNo == 1)  // "y"
        {
          assert(vectorIndex < values.size());
          values[vectorIndex++] = this->mesh_->meshPartition()->beginNodeGlobal(1) * this->meshWidth_ 
            + nodeY * this->meshWidth_;
        }
        else  // "z"
        {
          assert(vectorIndex < values.size());
          values[vectorIndex++] = this->mesh_->meshPartition()->beginNodeGlobal(2) * this->meshWidth_ 
            + nodeZ * this->meshWidth_;
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

//! for a specific component, get values from their local dof no.s
template<typename BasisOnMeshType, int nComponents>
template<int N>
void FieldVariableSetGetRegularFixed<BasisOnMeshType,nComponents>::
getValues(int componentNo, std::array<dof_no_t,N> dofLocalNo, std::array<double,N> &values)
{
  // if this is not a geometry field get the stored values
  if (!this->isGeometryField_)
  {
    FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
      template getValues<N>(componentNo, dofLocalNo, values);
    return;
  }

  // for geometry field compute information
  const node_no_t nLocalNodesInXDirection = this->mesh_->nNodesLocalWithGhosts(0);
  const node_no_t nLocalNodesInYDirection = this->mesh_->nNodesLocalWithGhosts(1);
  const int nDofsPerNode = BasisOnMeshType::nDofsPerNode();

  // loop over entries in values to be filled
  for (int i=0; i<N; i++)
  {
    int nodeLocalNo = int(dofLocalNo[i] / nDofsPerNode);
    int nodeLocalDofIndex = int(dofLocalNo[i] % nDofsPerNode);

    if (nodeLocalDofIndex > 0)   // if this is a derivative of Hermite, set to 0
    {
      values[i] = 0;
    }
    else
    {
      if (componentNo == 0)   // x direction
      {
        values[i] = this->mesh_->meshPartition()->beginNodeGlobal(0) * this->meshWidth_ 
        + (nodeLocalNo % nLocalNodesInXDirection) * this->meshWidth_;
      }
      else if (componentNo == 1)   // y direction
      {
        values[i] = this->mesh_->meshPartition()->beginNodeGlobal(1) * this->meshWidth_ 
        + (int(nodeLocalNo / nLocalNodesInXDirection) % nLocalNodesInYDirection) * this->meshWidth_;
      }
      else  // z direction
      {
        values[i] = this->mesh_->meshPartition()->beginNodeGlobal(2) * this->meshWidth_ 
        + int(nodeLocalNo / (nLocalNodesInXDirection*nLocalNodesInYDirection)) * this->meshWidth_;
      }
    }
  }
}
//! for a specific component, get values from their local dof no.s, as vector
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetRegularFixed<BasisOnMeshType,nComponents>::
getValues(int componentNo, std::vector<dof_no_t> dofLocalNo, std::vector<double> &values)
{
  // if this is not a geometry field get the stored values
  if (!this->isGeometryField_)
  {
    FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
      template getValues(componentNo, dofLocalNo, values);
    return;
  }

  // for geometry field compute information
  const node_no_t nLocalNodesInXDirection = this->mesh_->nNodesLocalWithGhosts(0);
  const node_no_t nLocalNodesInYDirection = this->mesh_->nNodesLocalWithGhosts(1);
  const int nDofsPerNode = BasisOnMeshType::nDofsPerNode();

  // resize result vector
  const int nValues = dofLocalNo.size();
  std::size_t previousSize = values.size();
  values.resize(previousSize+nValues);
  
  // loop over entries in values to be filled
  for (int i=0; i<nValues; i++)
  {
    int nodeLocalNo = int(dofLocalNo[i] / nDofsPerNode);
    int nodeLocalDofIndex = int(dofLocalNo[i] % nDofsPerNode);

    if (nodeLocalDofIndex > 0)   // if this is a derivative of Hermite, set to 0
    {
      values[previousSize+i] = 0;
    }
    else
    {
      if (componentNo == 0)   // x direction
      {
        values[previousSize+i] = this->mesh_->meshPartition()->beginNodeGlobal(0) * this->meshWidth_ 
        + (nodeLocalNo % nLocalNodesInXDirection) * this->meshWidth_;
      }
      else if (componentNo == 1)   // y direction
      {
        values[previousSize+i] = this->mesh_->meshPartition()->beginNodeGlobal(1) * this->meshWidth_ 
        + (int(nodeLocalNo / nLocalNodesInXDirection) % nLocalNodesInYDirection) * this->meshWidth_;
      }
      else  // z direction
      {
        values[previousSize+i] = this->mesh_->meshPartition()->beginNodeGlobal(2) * this->meshWidth_ 
        + int(nodeLocalNo / (nLocalNodesInXDirection*nLocalNodesInYDirection)) * this->meshWidth_;
      }
    }
  }
}

//! get values from their local dof no.s for all components
template<typename BasisOnMeshType, int nComponents>
template<int N>
void FieldVariableSetGetRegularFixed<BasisOnMeshType,nComponents>::
getValues(std::array<dof_no_t,N> dofLocalNo, std::array<std::array<double,nComponents>,N> &values)
{
  // if this is not a geometry field get the stored values
  if (!this->isGeometryField_)
  {
    FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
      template getValues<N>(dofLocalNo, values);
    return;
  }

  // for geometry field compute the entries
  const node_no_t nLocalNodesInXDirection = this->mesh_->nNodesLocalWithGhosts(0);
  const node_no_t nLocalNodesInYDirection = this->mesh_->nNodesLocalWithGhosts(1);
  const int nDofsPerNode = BasisOnMeshType::nDofsPerNode();

  // loop over entries in values to be filled
  for (int i=0; i<N; i++)
  {
    int nodeLocalNo = int(dofLocalNo[i] / nDofsPerNode);
    int nodeLocalDofIndex = int(dofLocalNo[i] % nDofsPerNode);

    if (nodeLocalDofIndex > 0)   // if this is a derivative of Hermite, set to 0
    {
      values[i][0] = 0;
      values[i][1] = 0;
      values[i][2] = 0;
    }
    else
    {
      // x direction
      values[i][0] = this->mesh_->meshPartition()->beginNodeGlobal(0) * this->meshWidth_ 
        + (nodeLocalNo % nLocalNodesInXDirection) * this->meshWidth_;

      // y direction
      values[i][1] = this->mesh_->meshPartition()->beginNodeGlobal(1) * this->meshWidth_ 
        + (int(nodeLocalNo / nLocalNodesInXDirection) % nLocalNodesInYDirection) * this->meshWidth_;

      // z direction
      values[i][2] = this->mesh_->meshPartition()->beginNodeGlobal(2) * this->meshWidth_ 
        + int(nodeLocalNo / (nLocalNodesInXDirection*nLocalNodesInYDirection)) * this->meshWidth_;
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
  std::array<dof_no_t,nDofsPerElement> elementDofs = this->mesh_->getElementDofLocalNos(elementNo);

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
  std::array<dof_no_t,nDofsPerElement> elementDofs = this->mesh_->getElementDofLocalNos(elementNo);

  // compute the corresponding geometry values
  this->getValues<nDofsPerElement>(elementDofs, values);
}

//! for a specific component, get a single value from local dof no.
template<typename BasisOnMeshType, int nComponents>
double FieldVariableSetGetRegularFixed<BasisOnMeshType,nComponents>::
getValue(int componentNo, node_no_t dofLocalNo)
{
  if (!this->isGeometryField_)
  {
    return FieldVariableSetGetStructured<BasisOnMeshType,nComponents>::
      getValue(componentNo, dofLocalNo);
  }

  // for geometry field compute information
  const node_no_t nLocalNodesInXDirection = this->mesh_->nNodesLocalWithGhosts(0);
  const node_no_t nLocalNodesInYDirection = this->mesh_->nNodesLocalWithGhosts(1);
  const int nDofsPerNode = BasisOnMeshType::nDofsPerNode();

  double value = 0;
  int nodeLocalNo = int(dofLocalNo / nDofsPerNode);
  int nodeLocalDofIndex = int(dofLocalNo % nDofsPerNode);

  if (nodeLocalDofIndex == 0)   // if this is not a derivative of Hermite (in which case it would be set to 0)
  {
    if (componentNo == 0)   // x direction
    {
      value = this->mesh_->meshPartition()->beginNodeGlobal(0) * this->meshWidth_ 
        + (nodeLocalNo % nLocalNodesInXDirection) * this->meshWidth_;
    }
    else if (componentNo == 1)   // y direction
    {
      value = this->mesh_->meshPartition()->beginNodeGlobal(1) * this->meshWidth_ 
        + (int(nodeLocalNo / nLocalNodesInXDirection) % nLocalNodesInYDirection) * this->meshWidth_;
    }
    else     // z direction
    {
      value = this->mesh_->meshPartition()->beginNodeGlobal(2) * this->meshWidth_ 
        + int(nodeLocalNo / (nLocalNodesInXDirection*nLocalNodesInYDirection)) * this->meshWidth_;
    }
  }
  return value;
}

//! copy the values from another field variable of the same type
template<typename BasisOnMeshType, int nComponents>
void FieldVariableSetGetRegularFixed<BasisOnMeshType,nComponents>::
setValues(FieldVariable<BasisOnMeshType,nComponents> &rhs)
{
  this->values_ = rhs.partitionedPetscVec();
}

};
