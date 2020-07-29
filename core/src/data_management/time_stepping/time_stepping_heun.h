#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "control/types.h"
#include "mesh/mesh.h"
#include "data_management/data.h"
#include "data_management/time_stepping/time_stepping.h"
#include "field_variable/field_variable.h"

namespace Data
{

/**  The datastructures used for Heun timestepping schemes, store all data as of timestepping
 *   schemes and an additional algebraic increment vector.
  */
template<typename FunctionSpaceType, int nComponents>
class TimeSteppingHeun : public TimeStepping<FunctionSpaceType,nComponents>
{
public:

  typedef FieldVariable::FieldVariable<FunctionSpaceType,nComponents> FieldVariableType;

  //! constructor
  TimeSteppingHeun(DihuContext context);

  //! destructur
  ~TimeSteppingHeun();

  //! return a reference to the algebraic increment vector, the PETSc Vec can be obtained via fieldVariable->valuesGlobal()
  std::shared_ptr<FieldVariableType> algebraicIncrement();

  //! print all stored data to stdout
  void print() override;

private:

  //! initializes the vectors with size
  void createPetscObjects() override;

  std::shared_ptr<FieldVariableType> algebraicIncrement_; //< the additional vector for delta u*
  // std::shared_ptr<FieldVariableType> algebraicSolution_; // / < the additional vector for u* // don't need this anymore
};

} // namespace Data

#include "data_management/time_stepping/time_stepping_heun.tpp"
