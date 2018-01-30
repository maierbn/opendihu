#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <array>
#include <map>

#include "field_variable/field_variable.h"
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

/** FieldVariable class for RegularFixed mesh
 */
template<int D, typename BasisFunctionType>
class FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>> :
  public FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>,
  public Interface<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>
{
public:
  typedef BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType> BasisOnMeshType;
 
  //! inherited constructor 
  using FieldVariableStructured<BasisOnMesh::BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>::FieldVariableStructured;
 
  //! set the meshWidths
  void setMeshWidth(std::array<double, D> &meshWidth);
  
  //! get the mesh width of a specific coordinate direction
  double meshWidth(int dimension) const;
  
  //! for a specific component, get values from their global dof no.s
  template<int N>
  void getValues(std::string component, std::array<dof_idx_t,N> dofGlobalNo, std::array<double,N> &values);
  
  //! get values from their global dof no.s for all components
  template<int N, int nComponents>
  void getValues(std::array<dof_idx_t,N> dofGlobalNo, std::array<std::array<double,nComponents>,N> &values);
    
  //! for a specific component, get the values corresponding to all element-local dofs
  template<int N>
  void getElementValues(std::string component, element_idx_t elementNo, std::array<double,BasisOnMeshType::nDofsPerElement()> &values);
  
  //! get the values corresponding to all element-local dofs for all components
  template<int nComponents>
  void getElementValues(element_idx_t elementNo, std::array<std::array<double,nComponents>,BasisOnMeshType::nDofsPerElement()> &values);
  
  //! for a specific component, get a single value from global dof no.
  double getValue(std::string component, node_idx_t dofGlobalNo);

  //! get a single value from global dof no. for all components
  template<int nComponents>
  std::array<double,nComponents> getValue(node_idx_t dofGlobalNo);

  //! write a exelem file header to a stream, for a particular element
  void outputHeaderExelem(std::ostream &file, element_idx_t currentElementGlobalNo, int fieldVariableNo=-1);

  //! write a exelem file header to a stream, for a particular element
  void outputHeaderExnode(std::ostream &file, node_idx_t currentNodeGlobalNo, int &valueIndex, int fieldVariableNo=-1);

  //! tell if 2 elements have the same exfile representation, i.e. same number of versions
  bool haveSameExfileRepresentation(element_idx_t element1, element_idx_t element2);

  //! get the internal PETSc vector values. The meaning of the values is instance-dependent (different for different BasisOnMeshTypes)
  Vec &values();
  
  //! get the number of components
  int nComponents() const;
  
  //! get the number of elements in the coordinate directions
  std::array<element_idx_t, BasisOnMeshType::Mesh::dim()> nElementsPerDimension() const;
  
  //! get the total number of elements
  element_idx_t nElements() const;
  
  //! get the names of the components that are part of this field variable
  std::vector<std::string> componentNames() const;
private:
 std::array<double, D> meshWidth_;   ///< the mesh width in each coordinate direction (has to be equal to work with FiniteElements)
};

};  // namespace

#include "field_variable/field_variable_regular_fixed.tpp"
