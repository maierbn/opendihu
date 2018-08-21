#include "field_variable/structured/04_field_variable_set_get_component_dependent_structured.h"

#include <sstream>
#include <cassert>
#include <array>
#include <type_traits>

namespace FieldVariable
{

//! get the values corresponding to all element-local dofs for all components
template<typename BasisOnMeshType>
void FieldVariableSetGetComponent<BasisOnMeshType,1>::
getElementValues(element_no_t elementNo, std::array<double,BasisOnMeshType::nDofsPerElement()> &values)
{
  VLOG(2) << "getElementValues element " << elementNo << ", 1 component";

  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();

  // prepare lookup indices for PETSc vector values_
  std::array<dof_no_t,nDofsPerElement> elementDofs = this->mesh_->getElementDofLocalNos(elementNo);

  this->values_->getValues(0, nDofsPerElement, (PetscInt *)elementDofs.data(), values.data());
}

//! get a single value from local dof no. for the single component
template<typename BasisOnMeshType>
double FieldVariableSetGetComponent<BasisOnMeshType,1>::
getValue(node_no_t dofLocalNo)
{
  double result;
  this->values_->getValues(0, 1, (PetscInt *)&dofLocalNo, &result);
  return result;
}


//! get a single value from local dof no. for all components
template<typename BasisOnMeshType, int nComponents>
std::array<double,nComponents> FieldVariableSetGetComponent<BasisOnMeshType,nComponents>::
getValue(node_no_t dofLocalNo)
{
  std::array<double,nComponents> result;

  // prepare lookup indices for PETSc vector values_
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    this->values_->getValues(componentNo, 1, &dofLocalNo, result.data()+componentNo);
  }

  return result;
}

//! get all stored local values, for 1 component
template<typename BasisOnMeshType>
void FieldVariableSetGetComponent<BasisOnMeshType,1>::
getValuesWithGhosts(std::vector<double> &values, bool onlyNodalValues)
{
  this->getValuesWithGhosts(0, values, onlyNodalValues);
}

//! get all stored local values, for 1 component
template<typename BasisOnMeshType>
void FieldVariableSetGetComponent<BasisOnMeshType,1>::
getValuesWithoutGhosts(std::vector<double> &values, bool onlyNodalValues)
{
  this->getValuesWithoutGhosts(0, values, onlyNodalValues);
}

//! set a single dof (all components), after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
template<typename BasisOnMeshType>
void FieldVariableSetGetComponent<BasisOnMeshType,1>::
setValue(dof_no_t dofLocalNo, double value, InsertMode petscInsertMode)
{
  this->values_->setValues(0, 1, (PetscInt*)&dofLocalNo, &value, petscInsertMode);
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. finishVectorManipulation must be called
}

//! set values for all components for dofs, after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
template<typename BasisOnMeshType>
void FieldVariableSetGetComponent<BasisOnMeshType,1>::
setValues(const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values, InsertMode petscInsertMode)
{
  this->values_->setValues(0, dofNosLocal.size(), (PetscInt*)dofNosLocal.data(), values.data(), petscInsertMode);
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. finishVectorManipulation must be called
}

//! set values for the single component for all local dofs, after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
template<typename BasisOnMeshType>
void FieldVariableSetGetComponent<BasisOnMeshType,1>::
setValuesWithGhosts(const std::vector<double> &values, InsertMode petscInsertMode)
{
  assert(values.size() == this->mesh_->meshPartition()->nDofsLocalWithGhosts());
  
  this->values_->setValues(0, values.size(), this->mesh_->meshPartition()->dofNosLocal().data(), values.data(), petscInsertMode);
}

//! set values for the single component for all local dofs, after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
template<typename BasisOnMeshType>
void FieldVariableSetGetComponent<BasisOnMeshType,1>::
setValuesWithoutGhosts(const std::vector<double> &values, InsertMode petscInsertMode)
{
  assert(values.size() == this->mesh_->meshPartition()->nDofsLocalWithoutGhosts());
  
  // set the values, this is the same call as setValuesWithGhosts, but the number of values is smaller and therefore the last dofs which are the ghosts are not touched
  this->values_->setValues(0, values.size(), this->mesh_->meshPartition()->dofNosLocal().data(), values.data(), petscInsertMode);
}

};  // namespace
