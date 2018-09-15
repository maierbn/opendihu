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

/**  The datastructures used for multidomain solver.
  */
template<typename FunctionSpaceType>
class Multidomain : public Data<FunctionSpaceType>
{
public:

  typedef FieldVariable::FieldVariable<FunctionSpaceType,3> GradientFieldVariableType;

  //! constructor of base class
  using Data<FunctionSpaceType,nComponents>::Data;

  //! return a reference to the rhs summand vector which is needed to apply the boundary conditions, the PETSc Vec can be obtained via fieldVariable->valuesGlobal()
  std::shared_ptr<GradientFieldVariableType> fibreDirection();

  //! print all stored data to stdout
  void print() override;

private:

  //! initializes the vectors with size
  void createPetscObjects() override;

  std::shared_ptr<FieldVariableType> fibreDirection_; ///< the direction of fibres
};

} // namespace Data

#include "data_management/multidomain.tpp"
