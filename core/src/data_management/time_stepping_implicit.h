#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "control/types.h"
#include "mesh/mesh.h"
#include "data_management/data.h"
#include "data_management/time_stepping.h"
#include "field_variable/field_variable.h"

namespace Data
{

/**  The datastructures used for implicit timestepping schemes.
  */
template<typename FunctionSpaceType, int nComponents>
class TimeSteppingImplicit : public TimeStepping<FunctionSpaceType,nComponents>
{
public:

  typedef FieldVariable::FieldVariable<FunctionSpaceType,nComponents> FieldVariableType;

  //! constructor of base class
  using TimeStepping<FunctionSpaceType,nComponents>::TimeStepping;

  //! return a reference to the rhs summand vector which is needed to apply the boundary conditions, the PETSc Vec can be obtained via fieldVariable->valuesGlobal()
  std::shared_ptr<FieldVariableType> boundaryConditionsRightHandSideSummand();

  //! return a reference to the rhs vector of the system, the PETSc Vec can be obtained via fieldVariable->valuesGlobal()
  std::shared_ptr<FieldVariableType> systemRightHandSide();

  //! print all stored data to stdout
  void print() override;

private:

  //! initializes the vectors with size
  void createPetscObjects() override;

  std::shared_ptr<FieldVariableType> boundaryConditionsRightHandSideSummand_; ///< what needs to be added to rhs when dirichlet boundary conditions are applied
  std::shared_ptr<FieldVariableType> systemRightHandSide_;                    ///< the rhs vector of the implicit system
};

} // namespace Data

#include "data_management/time_stepping_implicit.tpp"
