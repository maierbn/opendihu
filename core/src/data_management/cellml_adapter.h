#pragma once

#include <Python.h>  // has to be the first included header

#include "data_management/output_connector_data.h"
#include "control/python_config.h"

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
  typedef OutputConnectorData<FunctionSpaceType,nStates,nIntermediates> OutputConnectorDataType;

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

  //! give the names of all intermediates, will be called before initialize()
  void setIntermediateNames(const std::vector<std::string> &intermediateNames);

  //! return references to statesForTransfer_ and intermediatesForTransfer_
  void getStatesIntermediatesForTransfer(std::vector<int> &statesForTransfer, std::vector<int> &intermediatesForTransfer);

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the output_connector_data_transfer class
  std::shared_ptr<OutputConnectorDataType> getOutputConnectorData();

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

  //! initialize the outputConnectorData_ object, this can only be done after the states variable has been set
  void initializeOutputConnectorData();

  std::shared_ptr<FieldVariableIntermediates> intermediates_;   //< intermediates field variable
  std::shared_ptr<FieldVariableStates> states_;                 //< states field variable, this is a shared pointer with the timestepping scheme, which own the actual variable (creates it)
  std::vector<std::string> intermediateNames_;                  //< component names of the intermediates field variable

  std::shared_ptr<OutputConnectorDataType> outputConnectorData_;//< the object that holds all components of field variables that will be transferred to other solvers
  PythonConfig specificSettings_;                               //< the settings object

  std::vector<int> statesForTransfer_;                          //< state no.s to transfer to other solvers within output connector data
  std::vector<int> intermediatesForTransfer_;                   //< intermediate no.s to transfer to other solvers within output connector data

};

} // namespace Data

#include "data_management/cellml_adapter.tpp"
