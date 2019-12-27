#include "data_management/cellml_adapter.h"

namespace Data
{

template <int nStates, int nIntermediates, typename FunctionSpaceType>
CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::
CellmlAdapter(DihuContext context) :
  Data<FunctionSpaceType>::Data(context), specificSettings_(PythonConfig(context.getPythonConfig(), "CellML"))
{
}

template <int nStates, int nIntermediates, typename FunctionSpaceType>
void CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::
initialize()
{
  // call initialize of base class
  Data<FunctionSpaceType>::initialize();
}

template <int nStates, int nIntermediates, typename FunctionSpaceType>
void CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::
initializeOutputConnectorData()
{
  std::vector<int> statesForTransfer;
  std::vector<int> intermediatesForTransfer;

  this->specificSettings_.getOptionVector<int>("statesForTransfer", statesForTransfer);
  this->specificSettings_.getOptionVector<int>("intermediatesForTransfer", intermediatesForTransfer);

  LOG(DEBUG) << "parsed the following states for transfer: " << statesForTransfer << " (states: " << this->states() << "), intermediates: " << intermediatesForTransfer;

  // add states components
  for (std::vector<int>::iterator iter = statesForTransfer.begin(); iter != statesForTransfer.end(); iter++)
  {
    outputConnectorData_.addFieldVariable(this->states(), *iter);
  }

  // add intermediate components
  for (std::vector<int>::iterator iter = intermediatesForTransfer.begin(); iter != intermediatesForTransfer.end(); iter++)
  {
    outputConnectorData_.addFieldVariable2(this->intermediates(), *iter);
  }

  if (this->specificSettings_.hasKey("outputIntermediateIndex"))
  {
    LOG(WARNING) << specificSettings_ << "[\"outputIntermediateIndex\"] is no longer a valid option, use \"intermediatesForTransfer\" instead!";
  }

  if (this->specificSettings_.hasKey("outputStateIndex"))
  {
    LOG(WARNING) << specificSettings_ << "[\"outputStateIndex\"] is no longer a valid option, use \"statesForTransfer\" instead!";
  }
}

template <int nStates, int nIntermediates, typename FunctionSpaceType>
void CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::
setStatesVariable(std::shared_ptr<CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::FieldVariableStates> states)
{
  this->states_ = states;

  // after states variable has been set, continue initialize of output connector data
  initializeOutputConnectorData();
}

template <int nStates, int nIntermediates, typename FunctionSpaceType>
void CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::
createPetscObjects()
{
  this->intermediates_ = this->functionSpace_->template createFieldVariable<nIntermediates>("intermediates");
  this->intermediates_->setRepresentationContiguous();
}

template <int nStates, int nIntermediates, typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nIntermediates>> CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::
intermediates()
{
  return this->intermediates_;
}

template <int nStates, int nIntermediates, typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nStates>> CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::
states()
{
  return this->states_;
}

template <int nStates, int nIntermediates, typename FunctionSpaceType>
typename CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::OutputConnectorDataType &CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::
getOutputConnectorData()
{
  return this->outputConnectorData_;
}

template <int nStates, int nIntermediates, typename FunctionSpaceType>
typename CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::FieldVariablesForOutputWriter CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::
getFieldVariablesForOutputWriter()
{
  // these field variables will be written to output files
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> geometryField
    = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,3>>(this->functionSpace_->geometryField());

  return std::make_tuple(geometryField, intermediates_, states_);
}

} // namespace
