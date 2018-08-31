#include "field_variable/unstructured/04_field_variable_set_get_component_dependent_unstructured_deformable.h"

#include <sstream>
#include "utility/string_utility.h"
#include <map>
#include <fstream>
#include <iomanip>
#include <cassert>

namespace FieldVariable
{

using namespace StringUtility;


//! get a single value from local dof no. for all components
template<int D, typename BasisFunctionType, int nComponents>
std::array<double,nComponents> FieldVariableSetGet<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
getValue(node_no_t dofLocalNo) const
{
  std::array<double,nComponents> resultVector;

  for(int componentIndex = 0; componentIndex < nComponents; componentIndex++)
  {
    // get values and assign them to result values vector
    this->values_->getValues(componentIndex, 1, &dofLocalNo, resultVector.data() + componentIndex);
  }

  return resultVector;
}

//! get the values corresponding to all element-local dofs for all components
template<int D, typename BasisFunctionType>
void FieldVariableSetGet<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,1>::
getElementValues(element_no_t elementNo, std::array<double,FunctionSpaceType::nDofsPerElement()> &values) const
{
  assert(elementNo >= 0 && elementNo < this->functionSpace_->nElementsLocal());
  
  const int nDofsPerElement = FunctionSpaceType::nDofsPerElement();
  const std::vector<dof_no_t> &elementDofs = this->elementToDofMapping_->getElementDofs(elementNo);
  
  //VLOG(2) << "getElementValues element " << elementNo << ", nComponents=" << nComponents;

  this->values_->getValues(0, nDofsPerElement, (PetscInt *)elementDofs.data(), values.data());
}

//! get a single value from local dof no. for the single component
template<int D, typename BasisFunctionType>
double FieldVariableSetGet<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,1>::
getValue(node_no_t dofLocalNo) const
{
  double result;
  this->values_->getValues(0, 1, (PetscInt *)&dofLocalNo, &result);
  return result;
}

//! get all stored local values
template<int D, typename BasisFunctionType>
void FieldVariableSetGet<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,1>::
getValuesWithGhosts(std::vector<double> &values, bool onlyNodalValues) const
{
  this->getValuesWithGhosts(0, values, onlyNodalValues);
}

//! get all stored local values
template<int D, typename BasisFunctionType>
void FieldVariableSetGet<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,1>::
getValuesWithoutGhosts(std::vector<double> &values, bool onlyNodalValues) const
{
  this->getValuesWithoutGhosts(0, values, onlyNodalValues);
}

//! set a single dof (one components) , after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
template<int D, typename BasisFunctionType>
void FieldVariableSetGet<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,1>::
setValue(dof_no_t dofLocalNo, double value, InsertMode petscInsertMode)
{
  this->values_->setValues(0, 1, (PetscInt*)&dofLocalNo, &value, petscInsertMode);
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. finishVectorManipulation must be called
}

//! set values for one components for dofs, after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
template<int D, typename BasisFunctionType>
void FieldVariableSetGet<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,1>::
setValues(const std::vector<dof_no_t> &dofNosLocal, const std::vector<double> &values, InsertMode petscInsertMode)
{
  this->values_->setValues(0, dofNosLocal.size(), (PetscInt*)dofNosLocal.data(), values.data(), petscInsertMode);
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. finishVectorManipulation must be called
}

//! set values for the single component for all local dofs, after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
template<int D, typename BasisFunctionType>
void FieldVariableSetGet<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,1>::
setValuesWithGhosts(const std::vector<double> &values, InsertMode petscInsertMode)
{
  this->values_->setValues(0, values, petscInsertMode);
}

//! set values for the single component for all local dofs, after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
template<int D, typename BasisFunctionType>
void FieldVariableSetGet<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,1>::
setValuesWithoutGhosts(const std::vector<double> &values, InsertMode petscInsertMode)
{
  this->values_->setValues(0, values, petscInsertMode);
}

};  // namespace
