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


//! get a single value from global dof no. for all components
template<int D, typename BasisFunctionType, int nComponents>
std::array<double,nComponents> FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,nComponents>::
getValue(node_no_t dofGlobalNo)
{
  std::array<double,nComponents> resultVector;
  
  // transform global dof no.s to vector indices of first component
  PetscInt valuesVectorIndex = dofGlobalNo*nComponents;
  
  // create indices vector with values {0,1,2,...,nComponents-1}
  std::array<PetscInt,nComponents> indices;
  for(int i=0; i<nComponents; i++)
    indices[i] = valuesVectorIndex + i;
  
  // get values and assign them to result values vector
  VecGetValues(*this->values_, nComponents, indices.data(), resultVector.data());
  return resultVector;
}

//! get the values corresponding to all element-local dofs for all components
template<int D, typename BasisFunctionType>
void FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,1>::
getElementValues(element_no_t elementNo, std::array<double,BasisOnMeshType::nDofsPerElement()> &values)
{
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  const std::vector<dof_no_t> &elementDofs = this->elementToDofMapping_->getElementDofs(elementNo);
  VecGetValues(*this->values_, nDofsPerElement, (PetscInt *)elementDofs.data(), values.data());
}


//! get a single value from global dof no. for the single component
template<int D, typename BasisFunctionType>
double FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,1>::
getValue(node_no_t dofGlobalNo)
{
  double result;
  VecGetValues(*this->values_, 1, (PetscInt *)&dofGlobalNo, &result);
  return result;
}

//! set a single dof (all components) , after all calls to setValue(s), flushSetValues has to be called to apply the cached changes
template<int D, typename BasisFunctionType>
void FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,1>::
setValue(dof_no_t dofGlobalNo, double value, InsertMode petscInsertMode)
{
  VecSetValues(*this->values_, 1, (PetscInt*)&dofGlobalNo, &value, petscInsertMode);
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. flushSetValues must be called 
}

//! set values for all components for dofs, after all calls to setValue(s), flushSetValues has to be called to apply the cached changes
template<int D, typename BasisFunctionType>
void FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,1>::
setValues(std::vector<dof_no_t> &dofGlobalNos, std::vector<double> &values, InsertMode petscInsertMode)
{
  VecSetValues(*this->values_, dofGlobalNos.size(), (PetscInt*)dofGlobalNos.data(), values.data(), petscInsertMode);
  // after this VecAssemblyBegin() and VecAssemblyEnd(), i.e. flushSetValues must be called 
}

};  // namespace
