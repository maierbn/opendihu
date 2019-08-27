#pragma once

#include <Python.h>  // has to be the first included header

namespace Data
{

/**  The datastructures used for the cellml adapter.
  */
template <int nStates, int nIntermediates, typename FunctionSpaceType>
class CellmlAdapter : public Data<FunctionSpaceType>
{
public:
  typedef FieldVariable::FieldVariable<FunctionSpaceType,nIntermediates> FieldVariableIntermediates;
  typedef FieldVariable::FieldVariable<FunctionSpaceType,nStates> FieldVariableStates;

  //! constructor
  CellmlAdapter(DihuContext context);

  //! return a reference to the intermediates field variable
  std::shared_ptr<FieldVariableIntermediates> intermediates();

  //! return the field variable that stores all states
  std::shared_ptr<FieldVariableStates> states();

  //! initialize data structures
  void initialize();

  //! assign the field variable that stores the states, this variable is owned by the timestepping scheme
  void setStatesVariable(std::shared_ptr<FieldVariableStates> states);

  //! field variables that will be output by outputWriters
  typedef std::tuple<
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>>,     // geometry
    std::shared_ptr<FieldVariableIntermediates>,     // intermediates
    std::shared_ptr<FieldVariableStates>              // states
  > FieldVariablesForOutputWriter;

  //! get pointers to all field variables that can be written by output writers
  FieldVariablesForOutputWriter getFieldVariablesForOutputWriter();

private:

  //! initializes the vectors with size
  void createPetscObjects() override;

  std::shared_ptr<FieldVariableIntermediates> intermediates_;   ///< intermediates field variable
  std::shared_ptr<FieldVariableStates> states_;                 ///< states field variable, this is a shared pointer with the timestepping scheme, which own the actual variable (creates it)
};

} // namespace Data

#include "data_management/cellml_adapter.tpp"
