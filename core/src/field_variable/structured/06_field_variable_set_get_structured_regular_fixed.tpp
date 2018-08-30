#include "field_variable/structured/06_field_variable_set_get_structured_regular_fixed.h"

#include <sstream>
#include "easylogging++.h"
#include "utility/string_utility.h"

namespace FieldVariable
{

//! for a specific component, get all values
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetRegularFixed<FunctionSpaceType,nComponents>::
getValuesWithGhosts(int componentNo, std::vector<double> &values, bool onlyNodalValues)
{
  // if this is not a geometry field get the stored values
  if (!this->isGeometryField_)
  {
    VLOG(3) << "getValues(componentNo=" << componentNo << "): field variable is not a geometry field, retrieve stored values";
    FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
      getValuesWithGhosts(componentNo, values, onlyNodalValues);
    return;
  }

  // for geometry field compute information
  node_no_t nLocalNodesInXDirection = this->functionSpace_->nNodesLocalWithGhosts(0);
  node_no_t nLocalNodesInYDirection = 1;
  node_no_t nLocalNodesInZDirection = 1;

  const int D = FunctionSpaceType::dim();
  if (D >= 2)
    nLocalNodesInYDirection = this->functionSpace_->nNodesLocalWithGhosts(1);

  if (D >= 3)
    nLocalNodesInZDirection = this->functionSpace_->nNodesLocalWithGhosts(2);

  const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();

  // determine the number of values to be retrieved which is lower than the number of dofs for Hermite with only nodal values
  dof_no_t nValues = this->functionSpace_->nDofsLocalWithGhosts();
  if (onlyNodalValues)
  {
    nValues /= nDofsPerNode;
  }

  LOG(DEBUG) << "getValues, n dofs (with ghosts): " << this->functionSpace_->nDofsLocalWithGhosts() << ", nValues: " << nValues
    << ", nNodes: " << nLocalNodesInXDirection<< "," <<nLocalNodesInYDirection<< "," <<nLocalNodesInZDirection
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
          values[vectorIndex++] = this->functionSpace_->meshPartition()->beginNodeGlobal(0) * this->functionSpace_->meshWidth()
            + nodeX * this->functionSpace_->meshWidth();
        }
        else if (componentNo == 1)  // "y"
        {
          assert(vectorIndex < values.size());
          if (D < 2)
          {
            values[vectorIndex++] = 0;
          }
          else
          {
            values[vectorIndex++] = this->functionSpace_->meshPartition()->beginNodeGlobal(1) * this->functionSpace_->meshWidth()
              + nodeY * this->functionSpace_->meshWidth();
          }
        }
        else  // "z"
        {
          assert(vectorIndex < values.size());
          if (D < 3)
          {
            values[vectorIndex++] = 0;
          }
          else
          {
            values[vectorIndex++] = this->functionSpace_->meshPartition()->beginNodeGlobal(2) * this->functionSpace_->meshWidth()
              + nodeZ * this->functionSpace_->meshWidth();
          }
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

//! for a specific component, get all values
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetRegularFixed<FunctionSpaceType,nComponents>::
getValuesWithoutGhosts(int componentNo, std::vector<double> &values, bool onlyNodalValues)
{
  // if this is not a geometry field get the stored values
  if (!this->isGeometryField_)
  {
    VLOG(3) << "getValues(componentNo=" << componentNo << "): field variable is not a geometry field, retrieve stored values";
    FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
      getValuesWithoutGhosts(componentNo, values, onlyNodalValues);
    return;
  }

  // for geometry field compute information
  node_no_t nLocalNodesInXDirection = this->functionSpace_->nNodesLocalWithoutGhosts(0);
  node_no_t nLocalNodesInYDirection = 1;
  node_no_t nLocalNodesInZDirection = 1;

  const int D = FunctionSpaceType::dim();
  if (D >= 2)
    nLocalNodesInYDirection = this->functionSpace_->nNodesLocalWithoutGhosts(1);

  if (D >= 3)
    nLocalNodesInZDirection = this->functionSpace_->nNodesLocalWithoutGhosts(2);

  const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();

  // determine the number of values to be retrieved which is lower than the number of dofs for Hermite with only nodal values
  dof_no_t nValues = this->functionSpace_->nDofsLocalWithoutGhosts();
  if (onlyNodalValues)
  {
    nValues /= nDofsPerNode;
  }

  LOG(DEBUG) << "getValues, n dofs (without ghosts): " << this->functionSpace_->nDofsLocalWithoutGhosts() << ", nValues: " << nValues
    << ", nNodes: " << nLocalNodesInXDirection<< "," <<nLocalNodesInYDirection<< "," <<nLocalNodesInZDirection
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
          values[vectorIndex++] = this->functionSpace_->meshPartition()->beginNodeGlobal(0) * this->functionSpace_->meshWidth()
            + nodeX * this->functionSpace_->meshWidth();
        }
        else if (componentNo == 1)  // "y"
        {
          assert(vectorIndex < values.size());
          if (D < 2)
          {
            values[vectorIndex++] = 0;
          }
          else
          {
            values[vectorIndex++] = this->functionSpace_->meshPartition()->beginNodeGlobal(1) * this->functionSpace_->meshWidth()
              + nodeY * this->functionSpace_->meshWidth();
          }
        }
        else  // "z"
        {
          assert(vectorIndex < values.size());
          if (D < 3)
          {
            values[vectorIndex++] = 0;
          }
          else
          {
            values[vectorIndex++] = this->functionSpace_->meshPartition()->beginNodeGlobal(2) * this->functionSpace_->meshWidth()
              + nodeZ * this->functionSpace_->meshWidth();
          }
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
template<typename FunctionSpaceType, int nComponents>
template<int N>
void FieldVariableSetGetRegularFixed<FunctionSpaceType,nComponents>::
getValues(int componentNo, std::array<dof_no_t,N> dofLocalNo, std::array<double,N> &values)
{
  // if this is not a geometry field get the stored values
  if (!this->isGeometryField_)
  {
    FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
      template getValues<N>(componentNo, dofLocalNo, values);
    return;
  }
  
  // for geometry field compute information, this does not work for ghost dofs!
  node_no_t nLocalNodesInXDirection = this->functionSpace_->nNodesLocalWithoutGhosts(0);
  node_no_t nLocalNodesInYDirection = 1;

  const int D = FunctionSpaceType::dim();
  if (D >= 2)
    nLocalNodesInYDirection = this->functionSpace_->nNodesLocalWithoutGhosts(1);
  
  const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();

  // loop over entries in values to be filled
  for (int i=0; i<N; i++)
  {
    assert(dofLocalNo[i] < this->functionSpace_->nDofsLocalWithoutGhosts());
    
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
        values[i] = this->functionSpace_->meshPartition()->beginNodeGlobal(0) * this->functionSpace_->meshWidth()
        + (nodeLocalNo % nLocalNodesInXDirection) * this->functionSpace_->meshWidth();
      }
      else if (componentNo == 1)   // y direction
      {
        if (D < 2)
        {
          values[i] = 0;
        }
        else
        {
          values[i] = this->functionSpace_->meshPartition()->beginNodeGlobal(1) * this->functionSpace_->meshWidth()
            + (int(nodeLocalNo / nLocalNodesInXDirection) % nLocalNodesInYDirection) * this->functionSpace_->meshWidth();
        }
      }
      else  // z direction
      {
        if (D < 3)
        {
          values[i] = 0;
        }
        else
        {
          values[i] = this->functionSpace_->meshPartition()->beginNodeGlobal(2) * this->functionSpace_->meshWidth()
            + int(nodeLocalNo / (nLocalNodesInXDirection*nLocalNodesInYDirection)) * this->functionSpace_->meshWidth();
        }
      }
    }
  }
}
//! for a specific component, get values from their local dof no.s, as vector
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetRegularFixed<FunctionSpaceType,nComponents>::
getValues(int componentNo, std::vector<dof_no_t> dofLocalNo, std::vector<double> &values)
{
  // if this is not a geometry field get the stored values
  if (!this->isGeometryField_)
  {
    FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
      template getValues(componentNo, dofLocalNo, values);
    return;
  }

  // for geometry field compute information, this does not work for ghost dofs!
  node_no_t nLocalNodesInXDirection = this->functionSpace_->nNodesLocalWithoutGhosts(0);
  node_no_t nLocalNodesInYDirection = 1;

  const int D = FunctionSpaceType::dim();
  if (D >= 2)
    nLocalNodesInYDirection = this->functionSpace_->nNodesLocalWithoutGhosts(1);
  const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();

  // resize result vector
  const int nValues = dofLocalNo.size();
  std::size_t previousSize = values.size();
  values.resize(previousSize+nValues);
  
  // loop over entries in values to be filled
  for (int i=0; i<nValues; i++)
  {
    assert(dofLocalNo[i] < this->functionSpace_->nDofsLocalWithoutGhosts());
    
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
        values[previousSize+i] = this->functionSpace_->meshPartition()->beginNodeGlobal(0) * this->functionSpace_->meshWidth()
        + (nodeLocalNo % nLocalNodesInXDirection) * this->functionSpace_->meshWidth();
      }
      else if (componentNo == 1)   // y direction
      {
        if (D < 2)
        {
          values[previousSize+i] = 0;
        }
        else
        {
          values[previousSize+i] = this->functionSpace_->meshPartition()->beginNodeGlobal(1) * this->functionSpace_->meshWidth()
            + (int(nodeLocalNo / nLocalNodesInXDirection) % nLocalNodesInYDirection) * this->functionSpace_->meshWidth();
        }
      }
      else  // z direction
      {
        if (D < 3)
        {
          values[previousSize+i] = 0;
        }
        else
        {
          values[previousSize+i] = this->functionSpace_->meshPartition()->beginNodeGlobal(2) * this->functionSpace_->meshWidth()
            + int(nodeLocalNo / (nLocalNodesInXDirection*nLocalNodesInYDirection)) * this->functionSpace_->meshWidth();
        }
      }
    }
  }
}

//! get values from their local dof no.s for all components
template<typename FunctionSpaceType, int nComponents>
template<int N>
void FieldVariableSetGetRegularFixed<FunctionSpaceType,nComponents>::
getValues(std::array<dof_no_t,N> dofLocalNo, std::array<std::array<double,nComponents>,N> &values)
{
  // if this is not a geometry field get the stored values
  if (!this->isGeometryField_)
  {
    FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
      template getValues<N>(dofLocalNo, values);
    return;
  }

  // for geometry field compute the entries, this does not work for ghost dofs
  node_no_t nLocalNodesInXDirection = this->functionSpace_->nNodesLocalWithoutGhosts(0);
  node_no_t nLocalNodesInYDirection = 1;

  const int D = FunctionSpaceType::dim();
  if (D >= 2)
    nLocalNodesInYDirection = this->functionSpace_->nNodesLocalWithoutGhosts(1);
  const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();

  // loop over entries in values to be filled
  for (int i=0; i<N; i++)
  {
    assert(dofLocalNo[i] < this->functionSpace_->nDofsLocalWithoutGhosts());
    
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
      values[i][0] = this->functionSpace_->meshPartition()->beginNodeGlobal(0) * this->functionSpace_->meshWidth()
        + (nodeLocalNo % nLocalNodesInXDirection) * this->functionSpace_->meshWidth();

      // y direction
      if (D < 2)
      {
        values[i][1] = 0;
      }
      else
      {
        values[i][1] = this->functionSpace_->meshPartition()->beginNodeGlobal(1) * this->functionSpace_->meshWidth()
          + (int(nodeLocalNo / nLocalNodesInXDirection) % nLocalNodesInYDirection) * this->functionSpace_->meshWidth();
      }

      // z direction
      if (D < 3)
      {
        values[i][2] = 0;
      }
      else
      {
        values[i][2] = this->functionSpace_->meshPartition()->beginNodeGlobal(2) * this->functionSpace_->meshWidth()
          + int(nodeLocalNo / (nLocalNodesInXDirection*nLocalNodesInYDirection)) * this->functionSpace_->meshWidth();
      }
    }
  }
}

//! for a specific component, get the values corresponding to all element-local dofs
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetRegularFixed<FunctionSpaceType,nComponents>::
getElementValues(int componentNo, element_no_t elementNo, std::array<double,FunctionSpaceType::nDofsPerElement()> &values)
{
  // if this is not a geometry field get the stored values
  if (!this->isGeometryField_)
  {
    FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
      getElementValues(componentNo, elementNo, values);
    return;
  }

  // if this is a geometry field, compute the entries
  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();

  // get the element-local dofs of the element
  std::array<dof_no_t,nDofsPerElement> elementDofs = this->functionSpace_->getElementDofNosLocal(elementNo);

  // get the values
  this->getValues<nDofsPerElement>(componentNo, elementDofs, values);
}

//! get the values corresponding to all element-local dofs for all components
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetRegularFixed<FunctionSpaceType,nComponents>::
getElementValues(element_no_t elementNo, std::array<std::array<double,nComponents>,FunctionSpaceType::nDofsPerElement()> &values)
{
  // if this is not a geometry field get the stored values
  if (!this->isGeometryField_)
  {
    FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
      getElementValues(elementNo, values);
    return;
  }

  // for geometry field compute the requested values
  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();

  // get the element-local dofs of the element
  std::array<dof_no_t,nDofsPerElement> elementDofs = this->functionSpace_->getElementDofNosLocal(elementNo);

  // compute the corresponding geometry values
  this->getValues<nDofsPerElement>(elementDofs, values);
}

//! for a specific component, get a single value from local dof no.
template<typename FunctionSpaceType, int nComponents>
double FieldVariableSetGetRegularFixed<FunctionSpaceType,nComponents>::
getValue(int componentNo, node_no_t dofLocalNo)
{
  if (!this->isGeometryField_)
  {
    return FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
      getValue(componentNo, dofLocalNo);
  }

  assert(dofLocalNo < this->functionSpace_->meshPartition()->nDofsLocalWithoutGhosts());
    
  // for geometry field compute information, this does not work for ghost dofs
  node_no_t nLocalNodesInXDirection = this->functionSpace_->nNodesLocalWithoutGhosts(0);
  node_no_t nLocalNodesInYDirection = 1;

  const int D = FunctionSpaceType::dim();
  if (D >= 2)
    nLocalNodesInYDirection = this->functionSpace_->nNodesLocalWithoutGhosts(1);

  const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();

  double value = 0;
  int nodeLocalNo = int(dofLocalNo / nDofsPerNode);
  int nodeLocalDofIndex = int(dofLocalNo % nDofsPerNode);

  if (nodeLocalDofIndex == 0)   // if this is not a derivative of Hermite (in which case it would be set to 0)
  {
    if (componentNo == 0)   // x direction
    {
      value = this->functionSpace_->meshPartition()->beginNodeGlobal(0) * this->functionSpace_->meshWidth()
        + (nodeLocalNo % nLocalNodesInXDirection) * this->functionSpace_->meshWidth();
    }
    else if (componentNo == 1)   // y direction
    {
      if (D < 2)
      {
        value = 0;
      }
      else
      {
        value = this->functionSpace_->meshPartition()->beginNodeGlobal(1) * this->functionSpace_->meshWidth()
          + (int(nodeLocalNo / nLocalNodesInXDirection) % nLocalNodesInYDirection) * this->functionSpace_->meshWidth();
      }
    }
    else     // z direction
    {
      if (D < 3)
      {
        value = 0;
      }
      else
      {
        value = this->functionSpace_->meshPartition()->beginNodeGlobal(2) * this->functionSpace_->meshWidth()
          + int(nodeLocalNo / (nLocalNodesInXDirection*nLocalNodesInYDirection)) * this->functionSpace_->meshWidth();
      }
    }
  }
  return value;
}

//! copy the values from another field variable of the same type
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetRegularFixed<FunctionSpaceType,nComponents>::
setValues(FieldVariable<FunctionSpaceType,nComponents> &rhs)
{
  this->values_->setValues(*rhs.partitionedPetscVec());
}

};
