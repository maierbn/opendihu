#include "data_management/cellml_adapter.h"

namespace Data
{

template <int nStates, int nAlgebraics, typename FunctionSpaceType>
CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
CellmlAdapter(DihuContext context) :
  Data<FunctionSpaceType>::Data(context), specificSettings_(PythonConfig(context.getPythonConfig(), "CellML"))
{
}

//! constructor as a copy
template <int nStates, int nAlgebraics, typename FunctionSpaceType>
CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
CellmlAdapter(DihuContext context, const CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType> &rhs) :
  Data<FunctionSpaceType>::Data(context), specificSettings_(rhs.specificSettings_)
{
  // the followign variables are not copied:
  // std::shared_ptr<FieldVariableAlgebraics> algebraics_;   //< algebraics field variable
  // std::shared_ptr<FieldVariableStates> states_;           //< states field variable, this is a shared pointer with the timestepping scheme, which own the actual variable (creates it)
  // std::shared_ptr<FieldVariableAlgebraics> parameters_;   //< parameters field variable, the number of components is equal or less than the number of algebraics in order to not have to specify the number of parameters at compile time. This possibly creates a vector that is too large which is not harmful.
  // double *parameterValues_;                               //< a pointer to the data of the parameters_ Petsc Vec of the field variable

  algebraicNames_ = rhs.algebraicNames_;               //< component names of the algebraics field variable
  parameterNames_ = rhs.parameterNames_;               //< component names of the parameter field variable
  slotNames_ = rhs.slotNames_;                    //< names of the data slots that are used for slot connectors

  slotConnectorData_ = rhs.slotConnectorData_;    //< the object that holds all components of field variables that will be transferred to other solvers

  statesForTransfer_ = rhs.statesForTransfer_;                    //< state no.s to transfer to other solvers within slot connector data
  algebraicsForTransfer_ = rhs.algebraicsForTransfer_;                //< algebraic no.s to transfer to other solvers within slot connector data
  parametersForTransfer_ = rhs.parametersForTransfer_;                //< parameter no.s to transfer to other solvers within slot connector data
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
initializeSlotConnectorData()
{
  LOG(DEBUG) << "got the following states for transfer: " << statesForTransfer_ << " (states: " << this->states()
    << "), algebraics: " << algebraicsForTransfer_ << ", parameters: " << parametersForTransfer_;

  slotConnectorData_ = std::make_shared<SlotConnectorDataType>();

  // add states components
  for (std::vector<int>::iterator iter = statesForTransfer_.begin(); iter != statesForTransfer_.end(); iter++)
  {
    slotConnectorData_->addFieldVariable(this->states(), *iter);
  }

  // add algebraic components
  for (std::vector<int>::iterator iter = algebraicsForTransfer_.begin(); iter != algebraicsForTransfer_.end(); iter++)
  {
    slotConnectorData_->addFieldVariable2(this->algebraics(), *iter);
  }

  // add parameters components
  for (std::vector<int>::iterator iter = parametersForTransfer_.begin(); iter != parametersForTransfer_.end(); iter++)
  {
    slotConnectorData_->addFieldVariable2(this->parameters(), *iter);
  }

  // add slot names if given
  slotConnectorData_->slotNames.assign(slotNames_.begin(), slotNames_.end());

  // make sure that there are as many slot names as slots
  slotConnectorData_->slotNames.resize(slotConnectorData_->nSlots());

  LOG(DEBUG) << "set slot names: " << slotNames_ << " -> " << slotConnectorData_->slotNames;

  LOG(DEBUG) << "slotConnectorData: " << *slotConnectorData_;
}

template <int nStates, int nAlgebraics, typename FunctionSpaceType>
void CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
setStatesVariable(std::shared_ptr<CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::FieldVariableStates> states)
{
  // this will be called by the time stepping scheme after initialize()
  this->states_ = states;

  // after states variable has been set, continue initialize of slot connector data
  initializeSlotConnectorData();
}

template <int nStates, int nAlgebraics, typename FunctionSpaceType>
void CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
setAlgebraicAndParameterNames(const std::vector<std::string> &algebraicNames, const std::vector<std::string> &parameterNames, const std::vector<std::string> &slotNames)
{
  algebraicNames_ = algebraicNames;
  parameterNames_ = parameterNames;
  slotNames_ = slotNames;
}

template <int nStates, int nAlgebraics, typename FunctionSpaceType>
void CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
createPetscObjects()
{
  LOG(DEBUG) << "CellmlAdapter::createPetscObjects";

  // The states field variable is allocated by the timestepping class because this is the solution vector that the timestepping scheme operates on.
  // It gets then passed to this class by the call to setStatesVariable.
  // Therefore, here we only create the algebraics and the parameters field variables.
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
std::shared_ptr<typename CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::SlotConnectorDataType>
CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
getSlotConnectorData()
{
  return this->slotConnectorData_;
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
