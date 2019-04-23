#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <array>
#include <map>
#include <petscvec.h>
#include "easylogging++.h"

#include "field_variable/structured/02_field_variable_data_structured.h"
#include "partition/partitioned_petsc_vec/partitioned_petsc_vec.h"

namespace FieldVariable
{

/** Functionality to construct a 2D surface field variable from a 3D volume field variable, only for StructuredDeformableOfDimension<2>
 */
template<typename FunctionSpaceType, int nComponents_>
class FieldVariableDataStructuredForSurface :
  public FieldVariableDataStructured<FunctionSpaceType,nComponents_>
{
public:
  //! inherited constructor
  using FieldVariableDataStructured<FunctionSpaceType,nComponents_>::FieldVariableDataStructured;
};

template<typename BasisFunctionType, int nComponents>
class FieldVariableDataStructuredForSurface<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<2>,BasisFunctionType>,nComponents> :
  public FieldVariableDataStructured<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<2>,BasisFunctionType>,nComponents>
{
public:
  using FieldVariableDataStructured<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<2>,BasisFunctionType>,nComponents>::FieldVariableDataStructured;

  //! constructor from a 3D field variables.
  FieldVariableDataStructuredForSurface(FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunctionType>,nComponents> &rhs, Mesh::face_t face);
};

} // namespace

#include "field_variable/structured/02b_field_variable_data_structured_for_surface.tpp"
