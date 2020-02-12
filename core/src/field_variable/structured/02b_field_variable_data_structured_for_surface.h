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

  typedef FieldVariable<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunctionType>,nComponents> FieldVariable3D;

  //! constructor from a 3D field variable, this 2D field variable is at the surface given by face (only 2- or 2+ is possible). ownRankInvolvedInOutput is set to false if the own rank does not have any part of the data on this surface.
  FieldVariableDataStructuredForSurface(FieldVariable3D &rhs, Mesh::face_t face, bool &ownRankInvolvedInOutput);

  //! set values from 3D field variables
  void setValues(FieldVariable3D &rhs);

protected:

  // for a 3D numbering scheme with x*y*z numbers (e.g. dofs or ranks), get the surface dofs for face
  void getSurfaceNumbers(const std::array<node_no_t,3> size, int nDofsPerNode, Mesh::face_t face, std::vector<node_no_t> &surfaceNumbers);

  std::vector<dof_no_t> surfaceDofs_;   ///< local dof nos of the 3D field variable that specify the surface

};

} // namespace

#include "field_variable/structured/02b_field_variable_data_structured_for_surface.tpp"
