#include "field_variable/structured/04_field_variable_set_get_structured_deformable.h"

#include <sstream>
#include "utility/string_utility.h"
#include <cassert>

namespace FieldVariable
{
 
/*
//! for a specific component, get a single value from global dof no.
template<int D, typename BasisFunctionType>
int FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::
getValue(std::string component, node_no_t dofGlobalNo)
{
  return FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::getValue(component, dofGlobalNo);
}
*/
 
//! get a single value from global dof no. for all components
template<int D, typename BasisFunctionType>
template<std::size_t nComponents>
std::array<double,nComponents> FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::
getValue(node_no_t dofGlobalNo)
{
  // use the implementation of FieldVariableStructured
  return FieldVariableSetGetStructured<BasisOnMeshType>::template getValue<nComponents>(dofGlobalNo);
}

//! copy the values from another field variable of the same type
template<int D,typename BasisFunctionType>
void FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::
setValues(FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>> &rhs)
{
  VecCopy(rhs.values_, this->values_);
}

/*
//! set values for dofs
template<int D,typename BasisFunctionType>
template<std::size_t nComponents>
void FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::
setValues(std::vector<dof_no_t> &dofGlobalNos, std::vector<std::array<double,nComponents>> &values)
{
  if (!this->isGeometryField_)
  {
    FieldVariableSetGetStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::template setValues<nComponents>(dofGlobalNos, values);
  }
}*/

/*
//! set a single value
template<int D,typename BasisFunctionType>
template<std::size_t nComponents>
void FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::
setValue(dof_no_t dofGlobalNo, std::array<double,nComponents> &value)
{
  FieldVariableSetGetStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::template setValue<nComponents>(dofGlobalNo, value);
}
*/

//! calls PETSc functions to "assemble" the vector, i.e. flush the cached changes
template<int D,typename BasisFunctionType>
void FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::
flushSetValues()
{
  FieldVariableSetGetStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::flushSetValues();
}

};