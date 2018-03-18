#pragma once

#include <Python.h>  // has to be the first included header

#include "field_variable/structured/03_field_variable_data_structured_deformable.h"
#include "field_variable/field_variable_set_get.h"
#include "basis_on_mesh/05_basis_on_mesh.h"

namespace FieldVariable
{

/** FieldVariable class for StructuredDeformable mesh
 */
template<int D, typename BasisFunctionType>
class FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>> :
  public FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>
{
public:
  typedef BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType> BasisOnMeshType;
 
  //! inherited constructor 
  using FieldVariableData<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::FieldVariableData;
  
  //! for a specific component, get a single value from global dof no.
  double getValue(std::string component, node_no_t dofGlobalNo)
  {
    return FieldVariableSetGetStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::
      getValue(component, dofGlobalNo);
  }
  
  //! get a single value from global dof no. for all components
  template<std::size_t nComponents>
  std::array<double,nComponents> getValue(node_no_t dofGlobalNo);
  
  //! copy the values from another field variable of the same type
  void setValues(FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>> &rhs);
  
  //! set values for all components for dofs, after all calls to setValue(s), flushSetValues has to be called to apply the cached changes
  template<std::size_t nComponents>
  void setValues(std::vector<dof_no_t> &dofGlobalNos, std::vector<std::array<double,nComponents>> &values)
  {
    if (!this->isGeometryField_)
    {
      FieldVariableSetGetStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::
        template setValues<nComponents>(dofGlobalNos, values);
    }
  }  

  //! set a single dof (all components) , after all calls to setValue(s), flushSetValues has to be called to apply the cached changes
  template<std::size_t nComponents>
  void setValue(dof_no_t dofGlobalNo, std::array<double,nComponents> &value)
  {
    FieldVariableSetGetStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::
      template setValue<nComponents>(dofGlobalNo, value);
  }
  
  //! calls PETSc functions to "assemble" the vector, i.e. flush the cached changes
  void flushSetValues();
};

};  // namespace

#include "field_variable/structured/04_field_variable_set_get_structured_deformable.tpp"
