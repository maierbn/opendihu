#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include "utility/petsc_utility.h"
#include <memory>

#include "control/types.h"
#include "mesh/mesh.h"
#include "data_management/data.h"
#include "field_variable/field_variable.h"
#include "partition/partitioned_petsc_vec/partitioned_petsc_vec.h"
#include "partition/partitioned_petsc_mat/partitioned_petsc_mat.h"
#include "data_management/output_connector_data.h"

namespace Data
{

/**  The datastructures used for timestepping schemes, essentially stores a solution and an increment vector.
 *   nComponents is a constant that determines the number of components of each field variable.
 *   The time stepping works independent of that constant, it is only important for the output files, how many values there are associated with a single node.
 *   E.g. for a cellml model the nComponents should be set to the number of components of this model.
  */
template<typename FunctionSpaceType, int nComponents>
class TimeStepping :
  public Data<FunctionSpaceType>
{
public:

  typedef FieldVariable::FieldVariable<FunctionSpaceType,nComponents> FieldVariableType;
  typedef FieldVariable::FieldVariable<FunctionSpaceType,1> ScalarFieldVariableType;
  typedef OutputConnectorData<FunctionSpaceType,nComponents> OutputConnectorDataType;

  //! constructor
  TimeStepping(DihuContext context);

  //! destructor
  ~TimeStepping();

  //! return a reference to the solution vector, the PETSc Vec can be obtained via fieldVariable->valuesGlobal()
  std::shared_ptr<FieldVariableType> solution();

  //! return a reference to the increment vector, the PETSc Vec can be obtained via fieldVariable->valuesGlobal()
  std::shared_ptr<FieldVariableType> increment();
  
  //! set the names of the components for the solution field variable
  void setComponentNames(std::vector<std::string> componentNames);
  
  //! print all stored data to stdout
  virtual void print();

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the output_connector_data_transfer class
  std::shared_ptr<OutputConnectorDataType> getOutputConnectorData();

  //! field variables that will be output by outputWriters
  typedef std::tuple<
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>>,  // geometry
    std::shared_ptr<FieldVariableType>,  // solution
    std::vector<std::shared_ptr<ScalarFieldVariableType>>    // additional field variables that are not computed but transferred
  > FieldVariablesForOutputWriter;

  //! get pointers to all field variables that can be written by output writers
  FieldVariablesForOutputWriter getFieldVariablesForOutputWriter();

  //! output the given data for debugging
  std::string getString(std::shared_ptr<OutputConnectorDataType> data);

protected:

  //! initializes the vectors with size
  virtual void createPetscObjects();

  std::shared_ptr<FieldVariableType> solution_;     //< the vector of the variable of interest
  std::shared_ptr<FieldVariableType> increment_;    //< the vector for delta u, (note, this might be reduced in future to only a sub-part of the whole data vector if memory consumption is a problem)
  std::vector<std::string> componentNames_;         //< names of the components of the solution and increment variables
  std::vector<std::shared_ptr<ScalarFieldVariableType>> additionalFieldVariables_;   //< additional field variables that are not used for computation but can be passed from the discretizableInTime_ object to the enclosing solvers
  
  std::string debuggingName_;                       //< a name identifier only used for debugging

  std::shared_ptr<OutputConnectorDataType> outputConnectorData_;  //< the object that holds output connector data that will be transferred between solvers

private:
  //! get maximum number of expected non-zeros in the system matrix
  void getPetscMemoryParameters(int &nNonZerosDiagonal, int &nNonZerosOffdiagonal);

};

} // namespace Data

#include "data_management/time_stepping/time_stepping.tpp"
