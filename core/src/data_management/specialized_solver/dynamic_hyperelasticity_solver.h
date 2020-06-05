#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "data_management/data.h"
#include "field_variable/field_variable.h"

namespace Data
{

/**  The datastructures used for the dynamic solver of solid mechanics problems.
 *   It only contains field variables used for the timestepping (u,v,a),
 *   the static field variables are stored in a different Data class.
  */
template<typename FunctionSpaceType>
class DynamicHyperelasticitySolver :
  public Data<FunctionSpaceType>
{
public:

  typedef FieldVariable::FieldVariable<FunctionSpaceType,3> DisplacementsFieldVariableType;

  //! define the type of output connection variables, i.e. the values that will be transferred if the solver is part of a splitting or coupling scheme
  //! Two different field variables can be used: they must have the same function space but can have a different number of components. For example, for the CellMLAdapter, there are the "states" and the "algebraics" field variables.
  //! In this example, we use twice "1" as number of components, but you could, e.g. have OutputConnectorData<FunctionSpaceType,3,4>, etc.
  //! For each field variable you can transfer an abritrary subset of their components.
  typedef OutputConnectorData<FunctionSpaceType,1,3> OutputConnectorDataType;

  //! constructor
  DynamicHyperelasticitySolver(DihuContext context);

  //! field variable of u
  std::shared_ptr<DisplacementsFieldVariableType> displacements();

  //! field variable of v
  std::shared_ptr<DisplacementsFieldVariableType> velocities();

  //! ∂Wext
  std::shared_ptr<DisplacementsFieldVariableType> externalVirtualWorkDead();

  //! ∂Wint_stress
  std::shared_ptr<DisplacementsFieldVariableType> internalVirtualWork();

  //! ∂Wint_acceleration
  std::shared_ptr<DisplacementsFieldVariableType> accelerationTerm();

  //! field variables that will be output by outputWriters
  //type to use if there is no fiber direction field variable
  typedef std::tuple<
    std::shared_ptr<DisplacementsFieldVariableType>,  // current geometry field
    std::shared_ptr<DisplacementsFieldVariableType>,  // displacements
    std::shared_ptr<DisplacementsFieldVariableType>,  // velocities
    std::shared_ptr<DisplacementsFieldVariableType>,  // ∂Wint_stress
    std::shared_ptr<DisplacementsFieldVariableType>,  // ∂Wint_acceleration
    std::shared_ptr<DisplacementsFieldVariableType>   // ∂Wext
  >
  FieldVariablesForOutputWriter;

  //! get pointers to all field variables that can be written by output writers
  FieldVariablesForOutputWriter getFieldVariablesForOutputWriter();

protected:

  //! initializes the vectors with size
  void createPetscObjects() override;

  std::shared_ptr<DisplacementsFieldVariableType> displacements_;               //< the displacements, u
  std::shared_ptr<DisplacementsFieldVariableType> velocities_;                  //< v, the velocities
  std::shared_ptr<DisplacementsFieldVariableType> internalVirtualWork_;         //< ∂Wint_stress
  std::shared_ptr<DisplacementsFieldVariableType> accelerationTerm_;            //< ∂Wint_acceleration
  std::shared_ptr<DisplacementsFieldVariableType> externalVirtualWorkDead_;     //< ∂Wext_dead
};

} // namespace Data

#include "data_management/specialized_solver/dynamic_hyperelasticity_solver.tpp"
