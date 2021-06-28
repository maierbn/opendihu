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

/**  The datastructures used for "my new timestepping solver".
  */
template<typename FunctionSpaceType>
class DiffusionAdvectionSolver : public Data<FunctionSpaceType>
{
public:

  //! define the type of a scalar field variable (1 component), for convenience
  typedef FieldVariable::FieldVariable<FunctionSpaceType,1> ScalarFieldVariableType;

  //! define the type of a vector-valued field variable (3 component), for convenience
  typedef FieldVariable::FieldVariable<FunctionSpaceType,3> VectorFieldVariableType;

  //! define the type of output connection variables, i.e. the values that will be transferred if the solver is part of a splitting or coupling scheme
  //! Two different field variables can be used: they must have the same function space but can have a different number of components. For example, for the CellMLAdapter, there are the "states" and the "algebraics" field variables.
  //! In this example, we use twice "1" as number of components, but you could, e.g. have SlotConnectorData<FunctionSpaceType,3,4>, etc.
  //! For each field variable you can transfer an abritrary subset of their components.
  typedef SlotConnectorData<FunctionSpaceType,1,1> SlotConnectorDataType;

  //! constructor
  DiffusionAdvectionSolver(DihuContext context);

  //! return a reference to solution
  std::shared_ptr<ScalarFieldVariableType> solution();

  //! return a reference to increment
  std::shared_ptr<ScalarFieldVariableType> increment();

  //! return reference to the v matrix
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> vMatrix();

  //! get maximum number of expected non-zeros in stiffness matrix
  static void getPetscMemoryParameters(int &nNonZerosDiagonal, int &nNonZerosOffdiagonal);

  // ... define the field variables that you need


  //! initialize and create all variables
  void initialize();

  //! print all stored data to stdout
  void print();

  //! return the object that will be used to transfer values between solvers, in this case this includes only Vm
  std::shared_ptr<SlotConnectorDataType> getSlotConnectorData();

  //! field variables that will be output by outputWriters
  typedef std::tuple<
    std::shared_ptr<VectorFieldVariableType>,     // geometry, this always has to be the first field variable, such that the output writer knows the geometry of the mesh
    std::shared_ptr<ScalarFieldVariableType>,     // variableA,
    std::shared_ptr<ScalarFieldVariableType>      // variableB
    // ... add all field variables that you want to have in the output file
  > FieldVariablesForOutputWriter;

  //! get pointers to all field variables that can be written by output writers
  FieldVariablesForOutputWriter getFieldVariablesForOutputWriter();

private:

  //! create all field variables with their respective sizes, this will be called automatically within initialize by the base class
  void createPetscObjects() override;

  std::shared_ptr<ScalarFieldVariableType> solution_;   //< .. add a description of field variable A here!
  std::shared_ptr<ScalarFieldVariableType> increment_;   //< .. add a description of field variable B here!
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> vMatrix_;      //< the standard stiffness matrix of the finite element formulation, dofs with Dirichlet BCs will get the columns and rows cleared and diagonal set to 1.

  std::shared_ptr<SlotConnectorDataType> slotConnectorData_;    //< the object that stores all components of field variables that will be transferred to other solvers

  // define all needed field variables or other data
};

} // namespace Data

#include "data_management/specialized_solver/darcy/diffusion_advection_solver.tpp"
