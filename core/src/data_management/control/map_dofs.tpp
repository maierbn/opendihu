#include "data_management/control/multiple_instances.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <numeric>
#include <memory>

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "control/dihu_context.h"
#include "utility/petsc_utility.h"

namespace Data
{

template<typename FunctionSpaceType, typename NestedSolverType>
MapDofs<FunctionSpaceType, NestedSolverType>::
MapDofs(DihuContext context) : Data<FunctionSpaceType>(context)
{
}

//! create the additionalFieldVariables_
template<typename FunctionSpaceType, typename NestedSolverType>
void MapDofs<FunctionSpaceType, NestedSolverType>::
initialize(int nAdditionalFieldVariables, std::shared_ptr<typename NestedSolverType::OutputConnectorDataType> nestedSolverOutputConnectorData)
{
  nAdditionalFieldVariables_ = nAdditionalFieldVariables;

  // call initialize of base class
  Data<FunctionSpaceType>::initialize();

  // initialize output connector data
  outputConnectorData_ = std::make_shared<OutputConnectorDataType>();

  std::get<0>(*outputConnectorData_) = nestedSolverOutputConnectorData;
  std::get<1>(*outputConnectorData_) = std::make_shared<OutputConnectorData<FunctionSpaceType,1>>();

  // add all additional field variables
  for (int additionalFieldVariableNo = 0; additionalFieldVariableNo < nAdditionalFieldVariables_; additionalFieldVariableNo++)
  {
    std::get<1>(*outputConnectorData_)->addFieldVariable(additionalFieldVariables_[additionalFieldVariableNo]);
  }
}

template<typename FunctionSpaceType, typename NestedSolverType>
void MapDofs<FunctionSpaceType, NestedSolverType>::
createPetscObjects()
{
  for (int i = 0; i < nAdditionalFieldVariables_; i++)
  {
    std::stringstream name;
    name << "additionalFieldVariable" << i;
    std::shared_ptr<FieldVariableType> additionalFieldVariable = this->functionSpace_->template createFieldVariable<1>(name.str());

    additionalFieldVariables_.push_back(additionalFieldVariable);
  }
}

//! get a reference to the additional field variables
template<typename FunctionSpaceType, typename NestedSolverType>
std::vector<std::shared_ptr<typename MapDofs<FunctionSpaceType, NestedSolverType>::FieldVariableType>> &MapDofs<FunctionSpaceType, NestedSolverType>::additionalFieldVariables()
{
  return additionalFieldVariables_;
}

//! Get the data that will be transferred in the operator splitting or coupling to the other term of the splitting/coupling.
//! the transfer is done by the output_connector_data_transfer class
template<typename FunctionSpaceType, typename NestedSolverType>
std::shared_ptr<typename MapDofs<FunctionSpaceType, NestedSolverType>::OutputConnectorDataType> MapDofs<FunctionSpaceType, NestedSolverType>::
getOutputConnectorData()
{
  return outputConnectorData_;
}


template<typename FunctionSpaceType, typename NestedSolverType>
typename MapDofs<FunctionSpaceType, NestedSolverType>::FieldVariablesForOutputWriter MapDofs<FunctionSpaceType, NestedSolverType>::
getFieldVariablesForOutputWriter()
{
  return std::make_tuple(additionalFieldVariables_);
}

} // namespace
