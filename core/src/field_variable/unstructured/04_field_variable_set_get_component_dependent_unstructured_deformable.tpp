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
std::array<double,nComponents> FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
getValue(node_no_t dofLocalNo)
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
void FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,1>::
getElementValues(element_no_t elementNo, std::array<double,BasisOnMeshType::nDofsPerElement()> &values)
{
  assert(elementNo >= 0 && elementNo < this->mesh_->nLocalElements());
  
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  const std::vector<dof_no_t> &elementDofs = this->elementToDofMapping_->getElementDofs(elementNo);
  
  //VLOG(2) << "getElementValues element " << elementNo << ", nComponents=" << nComponents;

  this->values_->getValues(0, nDofsPerElement, (PetscInt *)elementDofs.data(), values.data());
}

//! get a single value from local dof no. for the single component
template<int D, typename BasisFunctionType>
double FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,1>::
getValue(node_no_t dofLocalNo)
{
  double result;
  this->values_->getValues(0, 1, (PetscInt *)&dofLocalNo, &result);
  return result;
}

//! get all stored local values
template<int D, typename BasisFunctionType>
void FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,1>::
getLocalValues(std::vector<double> &values)
{
  this->values_->getLocalValues(0, values);
}
  
//! set a single dof (one components) , after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
template<int D, typename BasisFunctionType>
void FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,1>::
setValue(dof_no_t dofLocalNo, double value, InsertMode petscInsertMode)
{
  this->values_->setValues(0, 1, (PetscInt*)&dofLocalNo, &value, petscInsertMode);
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. finishVectorManipulation must be called
}

//! set values for one components for dofs, after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
template<int D, typename BasisFunctionType>
void FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,1>::
setValues(std::vector<dof_no_t> &dofLocalNos, std::vector<double> &values, InsertMode petscInsertMode)
{
  this->values_->setValues(0, dofLocalNos.size(), (PetscInt*)dofLocalNos.data(), values.data(), petscInsertMode);
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. finishVectorManipulation must be called
}

//! set values for the single component for all local dofs, after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
template<int D, typename BasisFunctionType>
void FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,1>::
setValuesWithGhosts(std::vector<double> &values, InsertMode petscInsertMode)
{
  this->values_->setValues(0, values, petscInsertMode);
}

//! set values for the single component for all local dofs, after all calls to setValue(s), finishVectorManipulation has to be called to apply the cached changes
template<int D, typename BasisFunctionType>
void FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,1>::
setValuesWithoutGhosts(std::vector<double> &values, InsertMode petscInsertMode)
{
  this->values_->setValues(0, values, petscInsertMode);
}

};  // namespace
