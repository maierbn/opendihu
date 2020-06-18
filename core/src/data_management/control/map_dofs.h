#pragma once

#include <Python.h>  // has to be the first included header

namespace Data
{

/**  Additional field variables for the MapDofs class.
 */
template<typename FunctionSpaceType, typename NestedSolverType>
class MapDofs :
  public Data<FunctionSpaceType>
{
public:

  typedef FieldVariable::FieldVariable<FunctionSpaceType,1> FieldVariableType;

  //! Define the type of data that will be transferred between solvers when there is a coupling scheme.
  typedef std::tuple<
    std::shared_ptr<typename NestedSolverType::OutputConnectorDataType>,
    std::shared_ptr<OutputConnectorData<FunctionSpaceType,1>>
  > OutputConnectorDataType;

  //! constructor
  MapDofs(DihuContext context);

  //! create the additionalFieldVariables_
  void initialize(int nAdditionalFieldVariables, std::shared_ptr<typename NestedSolverType::OutputConnectorDataType> nestedSolverOutputConnectorData);

  //! field variables that will be output by outputWriters
  typedef std::tuple<std::vector<std::shared_ptr<FieldVariableType>>> FieldVariablesForOutputWriter;

  //! get a reference to the additional field variables
  std::vector<std::shared_ptr<FieldVariableType>> &additionalFieldVariables();

  //! get pointers to all field variables that can be written by output writers
  FieldVariablesForOutputWriter getFieldVariablesForOutputWriter();

  //! Get the data that will be transferred in the operator splitting or coupling to the other term of the splitting/coupling.
  //! the transfer is done by the output_connector_data_transfer class
  std::shared_ptr<OutputConnectorDataType> getOutputConnectorData();

protected:

  //! initializes the vectors with size
  virtual void createPetscObjects();

  int nAdditionalFieldVariables_;                                             //< number of additional field variables
  std::shared_ptr<OutputConnectorDataType> outputConnectorData_;              //< the object that holds all components of field variables that will be transferred to other solvers

  std::vector<std::shared_ptr<FieldVariableType>> additionalFieldVariables_;  //< the additional field variables that are created by MapDofs
};

} // namespace Data

#include "data_management/control/map_dofs.tpp"
