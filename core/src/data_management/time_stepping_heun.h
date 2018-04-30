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
 
/**  The datastructures used for Heun timestepping schemes, store all data as of timestepping
 *   schemes and an additional intermediate increment vector.
  */
template<typename BasisOnMeshType, int nComponents>
class TimeSteppingHeun : public TimeStepping<BasisOnMeshType,nComponents>
{
public:

  typedef FieldVariable::FieldVariable<BasisOnMeshType,nComponents> FieldVariableType;

  //! constructor
  TimeSteppingHeun(DihuContext context);
  
  //! destructur
  ~TimeSteppingHeun();

  //! return a reference to the intermediate increment vector, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariableType &intermediateIncrement();
  
  //! print all stored data to stdout
  void print() override;
    
private:
 
  //! initializes the vectors with size
  void createPetscObjects() override;
  
  std::shared_ptr<FieldVariableType> intermediateIncrement_; ///< the additional vector for delta u*,
  
};

} // namespace Data

#include "data_management/time_stepping_heun.tpp"
