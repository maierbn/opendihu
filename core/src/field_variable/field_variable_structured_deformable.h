#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <array>
#include <map>

#include "field_variable/field_variable_base.h"
#include "field_variable/field_variable_interface.h"
#include "field_variable/field_variable_structured.h"
#include "field_variable/component.h"
#include "basis_on_mesh/04_basis_on_mesh_nodes.h"
#include "mesh/unstructured_deformable.h"
#include "field_variable/element_to_node_mapping.h"
#include "field_variable/node_to_dof_mapping.h"

namespace FieldVariable
{

/** FieldVariable class for StructuredDeformable mesh
 */
template<int D, typename BasisFunctionType>
class FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>> :
  public FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>>
  //public Interface<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>
{
public:
  typedef BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType> BasisOnMeshType;
 
  //! inherited constructor 
  using FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>>::FieldVariableStructured;
  
  //! for a specific component, get a single value from global dof no.
  double getValue(std::string component, node_no_t dofGlobalNo)
  {
    return FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>>::getValue(component, dofGlobalNo);
  }
  
  //! get a single value from global dof no. for all components
  template<std::size_t nComponents>
  std::array<double,nComponents> getValue(node_no_t dofGlobalNo);
  
  //! copy the values from another field variable of the same type
  void setValues(FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>> &rhs);
  
  //! set values for all components for dofs, after all calls to setValue(s), flushSetValues has to be called to apply the cached changes
  template<std::size_t nComponents>
  void setValues(std::vector<dof_no_t> &dofGlobalNos, std::vector<std::array<double,nComponents>> &values)
  {
    if (!this->isGeometryField_)
    {
      FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>>::template setValues<nComponents>(dofGlobalNos, values);
    }
  }  

  //! set a single dof (all components) , after all calls to setValue(s), flushSetValues has to be called to apply the cached changes
  template<std::size_t nComponents>
  void setValue(dof_no_t dofGlobalNo, std::array<double,nComponents> &value)
  {
    FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>>::template setValue<nComponents>(dofGlobalNo, value);
  }
  
  //! calls PETSc functions to "assemble" the vector, i.e. flush the cached changes
  void flushSetValues();
  
  //! write a exelem file header to a stream, for a particular element
  void outputHeaderExelem(std::ostream &file, element_no_t currentElementGlobalNo);

  //! write a exelem file header to a stream, for a particular element
  void outputHeaderExnode(std::ostream &file, node_no_t currentNodeGlobalNo, int &valueIndex);

  //! tell if 2 elements have the same exfile representation, i.e. same number of versions
  bool haveSameExfileRepresentation(element_no_t element1, element_no_t element2);

  //! get the internal PETSc vector values. The meaning of the values is instance-dependent (different for different BasisOnMeshTypes)
  Vec &values();
  
  //! get the number of components
  int nComponents() const;
  
  //! get the names of the components
  std::vector<std::string> componentNames() const;
  
  //! get the number of elements in the coordinate directions
  std::array<element_no_t, BasisOnMeshType::Mesh::dim()> nElementsPerCoordinateDirection() const;
  
  //! get the total number of elements
  element_no_t nElements() const;
};

};  // namespace

#include "field_variable/field_variable_structured_deformable.tpp"
