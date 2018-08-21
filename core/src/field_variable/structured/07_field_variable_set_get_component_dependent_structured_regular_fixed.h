#pragma once

#include <Python.h>  // has to be the first included header

#include "field_variable/structured/06_field_variable_set_get_structured_regular_fixed.h"
#include "field_variable/field_variable_set_get.h"
#include "basis_on_mesh/basis_on_mesh.h"

namespace FieldVariable
{

/** FieldVariable class for RegularFixed mesh
 */
template<int D, typename BasisFunctionType, int nComponents>
class FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,nComponents> :
  public FieldVariableSetGetRegularFixed<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,nComponents>
{
public:
  typedef BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType> BasisOnMeshType;

  //! inherited constructor
  using FieldVariableSetGetRegularFixed<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,nComponents>::FieldVariableSetGetRegularFixed;
  using FieldVariableSetGetRegularFixed<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,nComponents>::getValue;

  //! get a single value from local dof no. for all components
  std::array<double,nComponents> getValue(node_no_t dofLocalNo);
};

/** Partial specialization for single component field variable
 */
template<int D, typename BasisFunctionType>
class FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,1> :
  public FieldVariableSetGetRegularFixed<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>,1>
{
public:
  typedef BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType> BasisOnMeshType;

  //! inherited constructor
  using FieldVariableSetGetRegularFixed<BasisOnMeshType,1>::FieldVariableSetGetRegularFixed;
  
};

};  // namespace

#include "field_variable/structured/07_field_variable_set_get_component_dependent_structured_regular_fixed.tpp"
