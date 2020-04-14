#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "control/types.h"
#include "mesh/mesh.h"
#include "data_management/data.h"
// #include "data_management/time_stepping/time_stepping.h"
#include "data_management/specialized_solver/multidomain.h"
#include "field_variable/field_variable.h"

namespace Data
{

/**  The datastructures used for "Parallel in TIme with Implicit Euler".
  */
template<typename FunctionSpaceType>
class PinTMD : public Data<FunctionSpaceType>
{
public:

  //! define the type of a scalar field variable (1 component), for convenience
  typedef FieldVariable::FieldVariable<FunctionSpaceType,1> ScalarFieldVariableType;

  //! define the type of a vector-valued field variable (3 component), for convenience
  typedef FieldVariable::FieldVariable<FunctionSpaceType,3> VectorFieldVariableType;

  //! define the type of output connection variables, i.e. the values that will be transferred if the solver is part of a splitting or coupling scheme
  //! Two different field variables can be used: they must have the same function space but can have a different number of components. For example, for the CellMLAdapter, there are the "states" and the "intermediates" field variables.
  //! In this example, we use twice "1" as number of components, but you could, e.g. have OutputConnectorData<FunctionSpaceType,3,4>, etc.
  //! For each field variable you can transfer an abritrary subset of their components.
  typedef OutputConnectorData<FunctionSpaceType,1,1> OutputConnectorDataType;

  //! constructor
  PinTMD(DihuContext context);

  //! return a reference to solution
  std::shared_ptr<ScalarFieldVariableType> solution();

  //! return a reference to fieldVariableB
  std::shared_ptr<ScalarFieldVariableType> fieldVariableB();

  // ... define the field variables that you need


  //! initialize and create all variables
  void initialize();

  //! initialize the solution field variable
  void setSolutionVariable(std::shared_ptr<ScalarFieldVariableType> solution);

  //! print all stored data to stdout
  void print();

  //! return the object that will be used to transfer values between solvers, in this case this includes only Vm
  std::shared_ptr<OutputConnectorDataType> getOutputConnectorData();

  //! field variables that will be output by outputWriters
  typedef std::tuple<
    std::shared_ptr<VectorFieldVariableType>,     // geometry, this always has to be the first field variable, such that the output writer knows the geometry of the mesh
    std::shared_ptr<ScalarFieldVariableType>,     // solution_,
    std::shared_ptr<ScalarFieldVariableType>      // variableB
    // ... add all field variables that you want to have in the output file
  > FieldVariablesForOutputWriter;

  //! get pointers to all field variables that can be written by output writers
  FieldVariablesForOutputWriter getFieldVariablesForOutputWriter();

private:

  //! create all field variables with their respective sizes, this will be called automatically within initialize by the base class
  void createPetscObjects() override;

  std::shared_ptr<ScalarFieldVariableType> solution_;         //< the solution of the laplace problem
  std::shared_ptr<ScalarFieldVariableType> fieldVariableB_;   //< .. add a description of field variable B here!

  std::shared_ptr<OutputConnectorDataType> outputConnectorData_;    ///< the object that stores all components of field variables that will be transferred to other solvers

  // define all needed field variables or other data
};

} // namespace Data

#include "data_management/specialized_solver/PinT_MD.tpp"
