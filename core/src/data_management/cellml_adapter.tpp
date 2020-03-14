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
  LOG(DEBUG) << "got the following states for transfer: " << statesForTransfer_ << " (states: " << this->states()
    << "), intermediates: " << intermediatesForTransfer_;

  outputConnectorData_ = std::make_shared<OutputConnectorDataType>();

  // add states components
  for (std::vector<int>::iterator iter = statesForTransfer_.begin(); iter != statesForTransfer_.end(); iter++)
  {
    outputConnectorData_->addFieldVariable(this->states(), *iter);
  }

  // add intermediate components
  for (std::vector<int>::iterator iter = intermediatesForTransfer_.begin(); iter != intermediatesForTransfer_.end(); iter++)
  {
    outputConnectorData_->addFieldVariable2(this->intermediates(), *iter);
  }

  // add parameters components
  for (std::vector<int>::iterator iter = parametersForTransfer_.begin(); iter != parametersForTransfer_.end(); iter++)
  {
    outputConnectorData_->addFieldVariable2(this->parameters(), *iter);
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
setIntermediateNames(const std::vector<std::string> &intermediateNames)
{
  intermediateNames_ = intermediateNames;
}

template <int nStates, int nIntermediates, typename FunctionSpaceType>
void CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::
createPetscObjects()
{
  // The states field variable is allocated by the timestepping class because this is the solution vector that the timestepping scheme operates on.
  // It gets then passed to this class by the call to setStatesVariable.
  // Therefore, here we only create the intermediates and the parameters field variables
  this->intermediates_ = this->functionSpace_->template createFieldVariable<nIntermediates>("intermediates", intermediateNames_);
  this->intermediates_->setRepresentationContiguous();

  std::vector<std::string> parameterNames(nIntermediates);
  for (int i = 0; i < nIntermediates; i++)
  {
    std::stringstream s;
    s << "parameter_" << i;
    parameterNames[i] = s.str();
  }

  this->parameters_ = this->functionSpace_->template createFieldVariable<nIntermediates>("parameters", parameterNames);
}

//! return a reference to the parameters vector
template <int nStates, int nIntermediates, typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nIntermediates>> CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::
parameters()
{
  return this->parameters_;
}

template <int nStates, int nIntermediates, typename FunctionSpaceType>
double *CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::
parameterValues()
{
  return this->parameterValues_;
}

//! get the parameteValues_ pointer from the parameters field variable, then the field variable can no longer be used until restoreParameterValues() gets called
template <int nStates, int nIntermediates, typename FunctionSpaceType>
void CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::
prepareParameterValues()
{
  this->parameters_->setRepresentationContiguous();
  PetscErrorCode ierr;
  ierr = VecGetArray(this->parameters_->getValuesContiguous(), &parameterValues_); CHKERRV(ierr);
}

//! restore the parameterValues_ pointer, such that the field variable can be used again
template <int nStates, int nIntermediates, typename FunctionSpaceType>
void CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::
restoreParameterValues()
{
  PetscErrorCode ierr;
  ierr = VecRestoreArray(this->parameters_->getValuesContiguous(), &parameterValues_); CHKERRV(ierr);
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
std::vector<int> &CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::
statesForTransfer()
{
  return statesForTransfer_;
}

template <int nStates, int nIntermediates, typename FunctionSpaceType>
std::vector<int> &CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::
intermediatesForTransfer()
{
  return intermediatesForTransfer_;
}

template <int nStates, int nIntermediates, typename FunctionSpaceType>
std::vector<int> &CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::
parametersForTransfer()
{
  return parametersForTransfer_;
}

template <int nStates, int nIntermediates, typename FunctionSpaceType>
std::shared_ptr<typename CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::OutputConnectorDataType>
CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::
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
