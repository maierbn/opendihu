#include "data_management/cellml_adapter.h"

namespace Data
{

template <int nStates, int nAlgebraics, typename FunctionSpaceType>
CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
CellmlAdapter(DihuContext context) :
  Data<FunctionSpaceType>::Data(context), specificSettings_(PythonConfig(context.getPythonConfig(), "CellML"))
{
}

template <int nStates, int nAlgebraics, typename FunctionSpaceType>
void CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
initialize()
{
  // call initialize of base class
  Data<FunctionSpaceType>::initialize();
}

template <int nStates, int nAlgebraics, typename FunctionSpaceType>
void CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
initializeOutputConnectorData()
{
  LOG(DEBUG) << "got the following states for transfer: " << statesForTransfer_ << " (states: " << this->states()
    << "), algebraics: " << algebraicsForTransfer_ << ", parameters: " << parametersForTransfer_;

  outputConnectorData_ = std::make_shared<OutputConnectorDataType>();

  // add states components
  for (std::vector<int>::iterator iter = statesForTransfer_.begin(); iter != statesForTransfer_.end(); iter++)
  {
    outputConnectorData_->addFieldVariable(this->states(), *iter);
  }

  // add algebraic components
  for (std::vector<int>::iterator iter = algebraicsForTransfer_.begin(); iter != algebraicsForTransfer_.end(); iter++)
  {
    outputConnectorData_->addFieldVariable2(this->algebraics(), *iter);
  }

  // add parameters components
  for (std::vector<int>::iterator iter = parametersForTransfer_.begin(); iter != parametersForTransfer_.end(); iter++)
  {
    outputConnectorData_->addFieldVariable2(this->parameters(), *iter);
  }
}

template <int nStates, int nAlgebraics, typename FunctionSpaceType>
void CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
setStatesVariable(std::shared_ptr<CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::FieldVariableStates> states)
{
  this->states_ = states;

  // after states variable has been set, continue initialize of output connector data
  initializeOutputConnectorData();
}

template <int nStates, int nAlgebraics, typename FunctionSpaceType>
void CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
setAlgebraicAndParameterNames(const std::vector<std::string> &algebraicNames, const std::vector<std::string> &parameterNames)
{
  algebraicNames_ = algebraicNames;
  parameterNames_ = parameterNames;
}

template <int nStates, int nAlgebraics, typename FunctionSpaceType>
void CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
createPetscObjects()
{
  // The states field variable is allocated by the timestepping class because this is the solution vector that the timestepping scheme operates on.
  // It gets then passed to this class by the call to setStatesVariable.
  // Therefore, here we only create the algebraics and the parameters field variables
  this->algebraics_ = this->functionSpace_->template createFieldVariable<nAlgebraics>("algebraics", algebraicNames_);
  this->algebraics_->setRepresentationContiguous();

  std::vector<std::string> parameterNames;

  for (int i = 0; i < parameterNames_.size(); i++)
  {
    std::stringstream s;
    s << "(P)" << parameterNames_[i];
    parameterNames.push_back(s.str());
  }

  for (int i = parameterNames.size(); i < nAlgebraics; i++)
  {
    std::stringstream s;
    s << "unusedParameter_" << i;
    parameterNames.push_back(s.str());
  }

  this->parameters_ = this->functionSpace_->template createFieldVariable<nAlgebraics>("parameters", parameterNames);
}

//! return a reference to the parameters vector
template <int nStates, int nAlgebraics, typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nAlgebraics>> CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
parameters()
{
  return this->parameters_;
}

template <int nStates, int nAlgebraics, typename FunctionSpaceType>
double *CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
parameterValues()
{
  return this->parameterValues_;
}

//! get the parameteValues_ pointer from the parameters field variable, then the field variable can no longer be used until restoreParameterValues() gets called
template <int nStates, int nAlgebraics, typename FunctionSpaceType>
void CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
prepareParameterValues()
{
  //LOG(DEBUG) << "parameters: " << *this->parameters_;

  this->parameters_->setRepresentationContiguous();
  PetscErrorCode ierr;
  Vec contiguousVec = this->parameters_->getValuesContiguous();
  ierr = VecGetArray(contiguousVec, &parameterValues_); CHKERRV(ierr);

#if 0
  PetscInt nValues;
  ierr = VecGetLocalSize(contiguousVec, &nValues); CHKERRV(ierr);
  LOG(DEBUG) << "parameter values has " << nValues << " entries: ";

  if (nValues == 11)
    counter++;
  if (counter == 22)
    LOG(FATAL) << "end, component 3 should not be -75 but -75.0059";
  //for (int i = 0; i < nValues; i++)
  //  LOG(DEBUG) << i << ": " << parameterValues_[i];
#endif
}

//! restore the parameterValues_ pointer, such that the field variable can be used again
template <int nStates, int nAlgebraics, typename FunctionSpaceType>
void CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
restoreParameterValues()
{
  PetscErrorCode ierr;
  ierr = VecRestoreArray(this->parameters_->getValuesContiguous(), &parameterValues_); CHKERRV(ierr);
}

template <int nStates, int nAlgebraics, typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nAlgebraics>> CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
algebraics()
{
  return this->algebraics_;
}

template <int nStates, int nAlgebraics, typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nStates>> CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
states()
{
  return this->states_;
}

template <int nStates, int nAlgebraics, typename FunctionSpaceType>
std::vector<int> &CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
statesForTransfer()
{
  return statesForTransfer_;
}

template <int nStates, int nAlgebraics, typename FunctionSpaceType>
std::vector<int> &CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
algebraicsForTransfer()
{
  return algebraicsForTransfer_;
}

template <int nStates, int nAlgebraics, typename FunctionSpaceType>
std::vector<int> &CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
parametersForTransfer()
{
  return parametersForTransfer_;
}

template <int nStates, int nAlgebraics, typename FunctionSpaceType>
std::shared_ptr<typename CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::OutputConnectorDataType>
CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
getOutputConnectorData()
{
  return this->outputConnectorData_;
}

template <int nStates, int nAlgebraics, typename FunctionSpaceType>
typename CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::FieldVariablesForOutputWriter CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
getFieldVariablesForOutputWriter()
{
  // these field variables will be written to output files
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> geometryField
    = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,3>>(this->functionSpace_->geometryField());

  return std::make_tuple(geometryField, algebraics_, states_, parameters_);
}

} // namespace
