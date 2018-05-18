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
  std::array<dof_no_t,nDofsPerElement> elementDofs = this->mesh_->getElementDofNos(elementNo);

  VecGetValues(this->values_, nDofsPerElement, (PetscInt *)elementDofs.data(), values.data());
}

//! get a single value from global dof no. for the single component
template<typename BasisOnMeshType>
double FieldVariableSetGetComponent<BasisOnMeshType,1>::
getValue(node_no_t dofGlobalNo)
{
  double result;
  VecGetValues(this->values_, 1, (PetscInt *)&dofGlobalNo, &result);
  return result;
}


//! get a single value from global dof no. for all components
template<typename BasisOnMeshType, int nComponents>
std::array<double,nComponents> FieldVariableSetGetComponent<BasisOnMeshType,nComponents>::
getValue(node_no_t dofGlobalNo)
{
  std::array<PetscInt,nComponents> indices;
  std::array<double,nComponents> result;

  // prepare lookup indices for PETSc vector values_
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    indices[componentNo] = dofGlobalNo*nComponents + componentNo;
  }

  VecGetValues(this->values_, nComponents, indices.data(), result.data());
  return result;
}

//! set a single dof (all components) , after all calls to setValue(s), flushSetValues has to be called to apply the cached changes
template<typename BasisOnMeshType>
void FieldVariableSetGetComponent<BasisOnMeshType,1>::
setValue(dof_no_t dofGlobalNo, double value, InsertMode petscInsertMode)
{
  VecSetValues(this->values_, 1, (PetscInt*)&dofGlobalNo, &value, petscInsertMode);
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. flushSetValues must be called
}

//! set values for all components for dofs, after all calls to setValue(s), flushSetValues has to be called to apply the cached changes
template<typename BasisOnMeshType>
void FieldVariableSetGetComponent<BasisOnMeshType,1>::
setValues(std::vector<dof_no_t> &dofGlobalNos, std::vector<double> &values, InsertMode petscInsertMode)
{
  VecSetValues(this->values_, dofGlobalNos.size(), (PetscInt*)dofGlobalNos.data(), values.data(), petscInsertMode);
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. flushSetValues must be called
}

};
