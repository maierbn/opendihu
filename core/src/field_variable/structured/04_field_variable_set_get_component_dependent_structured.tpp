#include "field_variable/structured/04_field_variable_set_get_component_dependent_structured.h"

#include <sstream>
#include <cassert>
#include <array>
#include <type_traits>

namespace FieldVariable
{

//! get the values corresponding to all element-local dofs for all components
template<typename FunctionSpaceType>
void FieldVariableSetGetComponent<FunctionSpaceType,1>::
getElementValues(element_no_t elementNoLocal, std::array<double,FunctionSpaceType::nDofsPerElement()> &values) const
{
  VLOG(2) << "getElementValues element " << elementNoLocal << ", 1 component";
  assert(this->values_);

  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();

  // prepare lookup indices for PETSc vector values_
  std::array<dof_no_t,nDofsPerElement> elementDofs = this->functionSpace_->getElementDofNosLocal(elementNoLocal);

  this->values_->getValues(0, nDofsPerElement, (PetscInt *)elementDofs.data(), values.data());
}

//! get the values corresponding to all element-local dofs for all components
template<typename FunctionSpaceType>
void FieldVariableSetGetComponent<FunctionSpaceType,1>::
getElementValues(Vc::int_v elementNoLocal, std::array<Vc::double_v,FunctionSpaceType::nDofsPerElement()> &values) const
{
  VLOG(2) << "getElementValues element " << elementNoLocal << ", 1 component";
  assert(this->values_);

  const int nVcComponents = Vc::double_v::size();
  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();

  std::array<PetscInt, nDofsPerElement*nVcComponents> indices;
  std::array<Vc::double_v, nDofsPerElement> result;

  for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
  {
    for (int vcComponent = 0; vcComponent < nVcComponents; vcComponent++)
    {
      if (elementNoLocal[vcComponent] == -1)
        indices[dofIndex*nVcComponents + vcComponent] = 0;     // set index to 0, then here the dof 0 is retrieved and also further used in computation but it is discarded later in setValue
      else
        indices[dofIndex*nVcComponents + vcComponent] = this->functionSpace_->getDofNo(elementNoLocal[vcComponent], dofIndex);
    }
  }

  // get the values for the current component
  this->values_->getValues(0, nDofsPerElement*nVcComponents, indices.data(), (double *)&result[nDofsPerElement]);

  // copy result to output values
  for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
  {
    values[dofIndex] = result[dofIndex];
  }
}

//! get a single value from local dof no. for the single component
template<typename FunctionSpaceType>
double FieldVariableSetGetComponent<FunctionSpaceType,1>::
getValue(node_no_t dofLocalNo) const
{
  assert(this->values_);
  double result;
  this->values_->getValues(0, 1, (PetscInt *)&dofLocalNo, &result);
  return result;
}


//! get a single value from local dof no. for all components
template<typename FunctionSpaceType, int nComponents>
std::array<double,nComponents> FieldVariableSetGetComponent<FunctionSpaceType,nComponents>::
getValue(node_no_t dofLocalNo) const
{
  assert(this->values_);
  std::array<double,nComponents> result;

  // prepare lookup indices for PETSc vector values_
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    this->values_->getValues(componentNo, 1, &dofLocalNo, result.data()+componentNo);
  }

  return result;
}

//! get all stored local values, for 1 component
template<typename FunctionSpaceType>
void FieldVariableSetGetComponent<FunctionSpaceType,1>::
getValues(const std::vector<dof_no_t> &dofLocalNo, std::vector<double> &values) const
{
  this->getValues(0, dofLocalNo, values);
}

//! get all stored local values, for 1 component
template<typename FunctionSpaceType>
void FieldVariableSetGetComponent<FunctionSpaceType,1>::
getValuesWithGhosts(std::vector<double> &values, bool onlyNodalValues) const
{
  this->getValuesWithGhosts(0, values, onlyNodalValues);
}

//! get all stored local values, for 1 component
template<typename FunctionSpaceType>
void FieldVariableSetGetComponent<FunctionSpaceType,1>::
getValuesWithoutGhosts(std::vector<double> &values, bool onlyNodalValues) const
{
  this->getValuesWithoutGhosts(0, values, onlyNodalValues);
}

//! set a single dof (all components), after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
template<typename FunctionSpaceType>
void FieldVariableSetGetComponent<FunctionSpaceType,1>::
setValue(dof_no_t dofLocalNo, double value, InsertMode petscInsertMode)
{
  assert(this->values_);
  assert(dofLocalNo < this->functionSpace_->meshPartition()->nDofsLocalWithGhosts());

  this->values_->setValues(0, 1, (PetscInt*)&dofLocalNo, &value, petscInsertMode);
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. finishGhostManipulation must be called
}

//! set a single dof (all components), after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
template<typename FunctionSpaceType>
void FieldVariableSetGetComponent<FunctionSpaceType,1>::
setValue(const Vc::int_v &dofLocalNo, const Vc::double_v &value, InsertMode petscInsertMode)
{
  assert(this->values_);
  /*
  std::array<double,Vc::double_v::size()> data;
  value.store(data.data());

  // store Vc vectors in order to get the raw memory
  std::array<int,Vc::double_v::size()> indices;
  dofLocalNo.store(indices.data());
*/
  // count number of non-negative indices in dofLocalNo, it is assumed that they occur all before the negative indices
  int nEntries = Vc::double_v::size() - Vc::isnegative(dofLocalNo).count();

  //this->values_->setValues(0, nEntries, indices.data(), data.data(), petscInsertMode);

  this->values_->setValues(0, nEntries, (PetscInt *)&dofLocalNo, (double *)&value, petscInsertMode);
}

//! set Vc::double_v::size() dofs with the same value
template<typename FunctionSpaceType>
void FieldVariableSetGetComponent<FunctionSpaceType,1>::
setValue(const Vc::int_v &dofLocalNo, double value, InsertMode petscInsertMode)
{
  assert(this->values_);
  assert(dofLocalNo[0] < this->functionSpace_->meshPartition()->nDofsLocalWithGhosts());

  std::array<double,Vc::double_v::size()> data;
  data.fill(value);
/*
  // store Vc vectors in order to get the raw memory
  std::array<int,Vc::double_v::size()> indices;
  dofLocalNo.store(indices.data());
*/
  // count number of non-negative indices in dofLocalNo, it is assumed that they occur all before the negative indices
  int nEntries = Vc::double_v::size() - Vc::isnegative(dofLocalNo).count();
/*
  Vc::double_v vectorizedValue = value;
  this->values_->setValues(0, nEntries, indices.data(), data.data(), petscInsertMode);
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. finishGhostManipulation must be called
  */

  this->values_->setValues(0, nEntries, (PetscInt *)&dofLocalNo, data.data(), petscInsertMode);
}

//! set the values from a petsc Vec
template<typename FunctionSpaceType>
void FieldVariableSetGetComponent<FunctionSpaceType,1>::
setValues(Vec petscVector)
{
  assert(this->values_);
  this->values_->setValues(0, petscVector);
}

//! set values for all components for dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
template<typename FunctionSpaceType>
void FieldVariableSetGetComponent<FunctionSpaceType,1>::
setValues(const std::vector<dof_no_t> &dofNosLocal, std::vector<double> &values, InsertMode petscInsertMode)
{
  assert(this->values_);
  this->values_->setValues(0, dofNosLocal.size(), (PetscInt*)dofNosLocal.data(), values.data(), petscInsertMode);
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. finishGhostManipulation must be called
}

//! set values for all components for dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
template<typename FunctionSpaceType>
template<int nValues>
void FieldVariableSetGetComponent<FunctionSpaceType,1>::
setValues(const std::array<dof_no_t,nValues> dofNosLocal, std::array<double,nValues> values, InsertMode petscInsertMode)
{
  assert(this->values_);
  this->values_->setValues(0, nValues, (PetscInt*)dofNosLocal.data(), values.data(), petscInsertMode);
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. finishGhostManipulation must be called
}

//! set values for the single component for all local dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
template<typename FunctionSpaceType>
void FieldVariableSetGetComponent<FunctionSpaceType,1>::
setValuesWithGhosts(const std::vector<double> &values, InsertMode petscInsertMode)
{
  if (values.size() != this->functionSpace_->meshPartition()->nDofsLocalWithGhosts())
  {
    LOG(ERROR) << "FieldVariable<1>::setValuesWithGhosts() values size: " << values.size() << ", nDofsWithGhosts: " << this->functionSpace_->meshPartition()->nDofsLocalWithGhosts();
  }
  assert(values.size() == this->functionSpace_->meshPartition()->nDofsLocalWithGhosts());
  assert(this->values_);
  
  this->values_->setValues(0, values.size(), this->functionSpace_->meshPartition()->dofNosLocal().data(), values.data(), petscInsertMode);
}

//! set values for the single component for all local dofs, after all calls to setValue(s), finishGhostManipulation has to be called to apply the cached changes
template<typename FunctionSpaceType>
void FieldVariableSetGetComponent<FunctionSpaceType,1>::
setValuesWithoutGhosts(const std::vector<double> &values, InsertMode petscInsertMode)
{
  assert(values.size() == this->functionSpace_->meshPartition()->nDofsLocalWithoutGhosts());
  assert(this->values_);
  
  // set the values, this is the same call as setValuesWithGhosts, but the number of values is smaller and therefore the last dofs which are the ghosts are not touched
  this->values_->setValues(0, values.size(), this->functionSpace_->meshPartition()->dofNosLocal().data(), values.data(), petscInsertMode);
}

} // namespace
