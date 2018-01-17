#pragma once

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
  double getValue(std::string component, node_idx_t dofGlobalNo)
  {
    return FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformable<D>,BasisFunctionType>>::getValue(component, dofGlobalNo);
  }
  
  //! get a single value from global dof no. for all components
  template<int nComponents>
  std::array<double,nComponents> getValue(node_idx_t dofGlobalNo);
  
  //! write a exelem file header to a stream, for a particular element
  void outputHeaderExelem(std::ostream &file, element_idx_t currentElementGlobalNo);

  //! write a exelem file header to a stream, for a particular element
  void outputHeaderExnode(std::ostream &file, node_idx_t currentNodeGlobalNo, int &valueIndex);

  //! tell if 2 elements have the same exfile representation, i.e. same number of versions
  bool haveSameExfileRepresentation(element_idx_t element1, element_idx_t element2);

  //! get the internal PETSc vector values. The meaning of the values is instance-dependent (different for different BasisOnMeshTypes)
  Vec &values();
  
  //! get the number of components
  int nComponents() const;
  
  //! get the number of elements
  std::array<int, BasisOnMeshType::Mesh::dim()> nElements() const;
};

};  // namespace

#include "field_variable/field_variable_structured_deformable.tpp"