#include "field_variable/structured/06_field_variable_set_get_structured_regular_fixed.h"

#include <sstream>
#include "easylogging++.h"
#include "utility/string_utility.h"

namespace FieldVariable
{

//! for a specific component, get all values
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetRegularFixed<FunctionSpaceType,nComponents>::
getValuesWithGhosts(int componentNo, std::vector<double> &values, bool onlyNodalValues) const
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
          values[vectorIndex++] = this->functionSpace_->meshPartition()->beginNodeGlobalNatural(0) * this->functionSpace_->meshWidth()
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
            values[vectorIndex++] = this->functionSpace_->meshPartition()->beginNodeGlobalNatural(1) * this->functionSpace_->meshWidth()
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
            values[vectorIndex++] = this->functionSpace_->meshPartition()->beginNodeGlobalNatural(2) * this->functionSpace_->meshWidth()
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
getValuesWithoutGhosts(int componentNo, std::vector<double> &values, bool onlyNodalValues) const
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
          values[vectorIndex++] = this->functionSpace_->meshPartition()->beginNodeGlobalNatural(0) * this->functionSpace_->meshWidth()
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
            values[vectorIndex++] = this->functionSpace_->meshPartition()->beginNodeGlobalNatural(1) * this->functionSpace_->meshWidth()
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
            values[vectorIndex++] = this->functionSpace_->meshPartition()->beginNodeGlobalNatural(2) * this->functionSpace_->meshWidth()
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

//! get all values
//! @param onlyNodalValues: if this is true, for Hermite only the non-derivative values are retrieved
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetRegularFixed<FunctionSpaceType,nComponents>::
getValuesWithGhosts(std::vector<std::array<double,nComponents>> &values, bool onlyNodalValues) const
{
  assert(this->values_);

  // determine the number of values to be retrived which is lower than the number of dofs for Hermite with only nodal values
  dof_no_t nValues = this->functionSpace_->meshPartition()->nDofsLocalWithGhosts();
  if (onlyNodalValues)
  {
    const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();
    nValues /= nDofsPerNode;
  }

  // resize output vector
  VLOG(2) << "Field variable structured, getValues, resize values vector to " << nValues << " entries.";

  std::size_t previousSize = values.size();
  values.resize(previousSize+nValues);

  // loop over components and get data component-wise
  std::vector<double> buffer;
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    // get values into buffer
    getValuesWithGhosts(componentNo, buffer, onlyNodalValues);
    assert(nValues == buffer.size());

    // copy values from buffer to output vector
    for (int valueIndex = 0; valueIndex < nValues; valueIndex++)
    {
      values[previousSize+valueIndex][componentNo] = buffer[valueIndex];
    }
  }
}

//! get all values
//! @param onlyNodalValues: if this is true, for Hermite only the non-derivative values are retrieved
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetRegularFixed<FunctionSpaceType,nComponents>::
getValuesWithoutGhosts(std::vector<std::array<double,nComponents>> &values, bool onlyNodalValues) const
{
  // determine the number of values to be retrived which is lower than the number of dofs for Hermite with only nodal values
  dof_no_t nValues = this->functionSpace_->meshPartition()->nDofsLocalWithoutGhosts();
  if (onlyNodalValues)
  {
    const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();
    nValues /= nDofsPerNode;
  }

  // resize output vector
  std::size_t previousSize = values.size();
  values.resize(previousSize+nValues);
  VLOG(2) << "Field variable regular fixed, getValues, resize values vector to " << previousSize+nValues << " entries.";

  // loop over components and get data component-wise
  std::vector<double> buffer;
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    // get values into buffer
    getValuesWithoutGhosts(componentNo, buffer, onlyNodalValues);
    assert(nValues == buffer.size());

    // copy values from buffer to output vector
    for (int valueIndex = 0; valueIndex < nValues; valueIndex++)
    {
      values[previousSize+valueIndex][componentNo] = buffer[valueIndex];
    }
  }
}

//! for a specific component, get values from their local dof no.s
template<typename FunctionSpaceType, int nComponents>
template<int N>
void FieldVariableSetGetRegularFixed<FunctionSpaceType,nComponents>::
getValues(int componentNo, std::array<dof_no_t,N> dofLocalNo, std::array<double,N> &values) const
{
  // if this is not a geometry field get the stored values
  if (!this->isGeometryField_)
  {
    FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
      template getValues<N>(componentNo, dofLocalNo, values);
    return;
  }
  
  // for geometry field compute information
  const int D = FunctionSpaceType::dim();
  const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();

  // loop over entries in values to be filled
  for (int i=0; i<N; i++)
  {
    int nodeNoLocal = int(dofLocalNo[i] / nDofsPerNode);
    int nodeLocalDofIndex = int(dofLocalNo[i] % nDofsPerNode);
    std::array<global_no_t,D> coordinates = this->functionSpace_->meshPartition()->getCoordinatesGlobal(nodeNoLocal);

    if (nodeLocalDofIndex > 0)   // if this is a derivative of Hermite, set to 0
    {
      values[i] = 0;
    }
    else
    {
      if (componentNo == 0)   // x direction
      {
        values[i] = coordinates[0] * this->functionSpace_->meshWidth();;
      }
      else if (componentNo == 1)   // y direction
      {
        if (D < 2)
        {
          values[i] = 0;
        }
        else
        {
          values[i] = coordinates[1] * this->functionSpace_->meshWidth();;
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
          values[i] = coordinates[2] * this->functionSpace_->meshWidth();;
        }
      }
    }
  }
}
//! for a specific component, get values from their local dof no.s, as vector
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetRegularFixed<FunctionSpaceType,nComponents>::
getValues(int componentNo, const std::vector<dof_no_t> &dofLocalNo, std::vector<double> &values) const
{
  // if this is not a geometry field get the stored values
  if (!this->isGeometryField_)
  {
    FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::getValues(componentNo, dofLocalNo, values);
    return;
  }

  // for geometry field compute information
  const int D = FunctionSpaceType::dim();
  const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();

  // resize result vector
  const int nValues = dofLocalNo.size();
  std::size_t previousSize = values.size();
  values.resize(previousSize+nValues);
  
  // loop over entries in values to be filled
  for (int i=0; i<nValues; i++)
  {
    int nodeNoLocal = int(dofLocalNo[i] / nDofsPerNode);
    int nodeLocalDofIndex = int(dofLocalNo[i] % nDofsPerNode);
    std::array<global_no_t,D> coordinates = this->functionSpace_->meshPartition()->getCoordinatesGlobal(nodeNoLocal);

    if (nodeLocalDofIndex > 0)   // if this is a derivative of Hermite, set to 0
    {
      values[previousSize+i] = 0;
    }
    else
    {
      if (componentNo == 0)   // x direction
      {
        values[previousSize+i] = coordinates[0] * this->functionSpace_->meshWidth();
      }
      else if (componentNo == 1)   // y direction
      {
        if (D < 2)
        {
          values[previousSize+i] = 0;
        }
        else
        {
          values[previousSize+i] = coordinates[1] * this->functionSpace_->meshWidth();
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
          values[previousSize+i] = coordinates[2] * this->functionSpace_->meshWidth();
        }
      }
    }
  }
}

//! get values from their local dof no.s for all components
template<typename FunctionSpaceType, int nComponents>
template<int N>
void FieldVariableSetGetRegularFixed<FunctionSpaceType,nComponents>::
getValues(std::array<dof_no_t,N> dofLocalNo, std::array<std::array<double,nComponents>,N> &values) const
{
  // if this is not a geometry field get the stored values
  if (!this->isGeometryField_)
  {
    FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
      template getValues<N>(dofLocalNo, values);
    return;
  }

  // for geometry field compute the entries
  const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();
  const int D = FunctionSpaceType::dim();

  // loop over entries in values to be filled
  for (int i=0; i<N; i++)
  {
    int nodeNoLocal = int(dofLocalNo[i] / nDofsPerNode);
    int nodeLocalDofIndex = int(dofLocalNo[i] % nDofsPerNode);
    std::array<global_no_t,D> coordinates = this->functionSpace_->meshPartition()->getCoordinatesGlobal(nodeNoLocal);

    if (nodeLocalDofIndex > 0)   // if this is a derivative of Hermite, set to 0
    {
      values[i][0] = 0;
      values[i][1] = 0;
      values[i][2] = 0;
    }
    else
    {
      // x direction
      values[i][0] = coordinates[0] * this->functionSpace_->meshWidth();

      // y direction
      if (D < 2)
      {
        values[i][1] = 0;
      }
      else
      {
        values[i][1] = coordinates[1] * this->functionSpace_->meshWidth();
      }

      // z direction
      if (D < 3)
      {
        values[i][2] = 0;
      }
      else
      {
        values[i][2] = coordinates[2] * this->functionSpace_->meshWidth();
      }
    }
  }
}

//! for a specific component, get the values corresponding to all element-local dofs
template<typename FunctionSpaceType, int nComponents>
void FieldVariableSetGetRegularFixed<FunctionSpaceType,nComponents>::
getElementValues(int componentNo, element_no_t elementNo, std::array<double,FunctionSpaceType::nDofsPerElement()> &values) const
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
getElementValues(element_no_t elementNo, std::array<std::array<double,nComponents>,FunctionSpaceType::nDofsPerElement()> &values) const
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
getValue(int componentNo, node_no_t dofLocalNo) const
{
  if (!this->isGeometryField_)
  {
    return FieldVariableSetGetStructured<FunctionSpaceType,nComponents>::
      getValue(componentNo, dofLocalNo);
  }
    
  // for geometry field compute information
  const int D = FunctionSpaceType::dim();
  const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();
  int nodeNoLocal = int(dofLocalNo / nDofsPerNode);
  int nodeLocalDofIndex = int(dofLocalNo % nDofsPerNode);

  std::array<global_no_t,D> coordinates = this->functionSpace_->meshPartition()->getCoordinatesGlobal(nodeNoLocal);

  double value = 0;
  if (nodeLocalDofIndex == 0)   // if this is not a derivative of Hermite (in which case it would be set to 0)
  {
    if (componentNo == 0)   // x direction
    {
      value = coordinates[0] * this->functionSpace_->meshWidth();
    }
    else if (componentNo == 1)   // y direction
    {
      if (D < 2)
      {
        value = 0;
      }
      else
      {
        value = coordinates[1] * this->functionSpace_->meshWidth();
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
        value = coordinates[2] * this->functionSpace_->meshWidth();
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

}  // namespace
