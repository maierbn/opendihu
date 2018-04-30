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
 
  // sollte mitgeerbt werden:
  // typedef FieldVariable::FieldVariable<BasisOnMeshType,nComponents> FieldVariableType;
 
  //! constructor
  TimeSteppingHeun(DihuContext context);
  
  //! destructur
  ~TimeSteppingHeun();
  
  // return a reference to the solution vector, the PETSc Vec can be obtained via fieldVariable.values()
  // FieldVariableType &solution();
  
  // ! return a reference to the increment vector, the PETSc Vec can be obtained via fieldVariable.values()
  // FieldVariableType &increment();
  
  //! return a reference to the intermediate increment vector, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariableType &intermediate_increment();
  
  //! print all stored data to stdout
  void print() override;
  
  // ! return the number of degrees of freedom per mesh node
  //static constexpr int getNDofsPerNode();
  
  // ! return the total number of degrees of freedom, this can be a multiple of the number of nodes of the mesh
  //virtual dof_no_t nUnknowns();
  
  // ! field variables that will be output by outputWriters
  // typedef std::tuple<
  //  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>>,  // geometry
  //  std::shared_ptr<FieldVariableType>  // solution
  // > OutputFieldVariables;
  
  // ! get pointers to all field variables that can be written by output writers
  // OutputFieldVariables getOutputFieldVariables();
  
private:
 
  //! initializes the vectors with size
  void createPetscObjects() override;
  
  std::shared_ptr<FieldVariableType> intermediate_increment_; ///< the additional vector for delta u*,
  
};

} // namespace Data

#include "data_management/time_stepping_heun.tpp"
