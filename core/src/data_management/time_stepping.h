#pragma once

#include <petscmat.h>
#include <memory>

#include "control/types.h"
#include "mesh/mesh.h"
#include "data_management/data.h"
#include "field_variable/field_variable.h"

class DihuContext;

namespace Data
{
 
template<typename BasisOnMeshType>
class TimeStepping : public Data<BasisOnMeshType>
{
public:
 
  //! constructor
  TimeStepping(const DihuContext &context);
  
  //! destructur
  ~TimeStepping();
  
  //! return a reference to the solution vector, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<BasisOnMeshType> &solution();
  
  //! return a reference to the increment vector, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<BasisOnMeshType> &increment();
 
  //! print all stored data to stdout
  void print();
  
private:
 
  //! initializes the vectors with size
  void createPetscObjects();
  
  bool disablePrinting_ = false;    ///< if printing vectors is disabled
  
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>> solution_;            ///< the vector of the variable of interest
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>> increment_;        ///< the vector for delta u
};

} // namespace Data

#include "data_management/time_stepping.tpp"