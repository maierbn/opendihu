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

  slotConnectorDataTimestepping_ = rhs.slotConnectorDataTimestepping_;    //< the object that holds all components of field variables that will be transferred to other solvers

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
setStatesVariable(std::shared_ptr<CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::FieldVariableStates> states)
{
  // this will be called by the time stepping scheme after initialize() and before setSlotConnectorDataTimestepping()
  this->states_ = states;
}

template <int nStates, int nAlgebraics, typename FunctionSpaceType>
void CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
setSlotConnectorDataTimestepping(std::shared_ptr<SlotConnectorData<FunctionSpaceType,nStates>> slotConnectorDataTimestepping)
{
  slotConnectorDataTimestepping_ = slotConnectorDataTimestepping;

  // This method is called once in initialize() of the timestepping scheme.
  // Add all state and algebraic values for transfer (option "algebraicsForTransfer"), which are stored in this->data_.getSlotConnectorData().
  // The states are store in variable1 of slotConnectorDataTimestepping, the algebraics are stored in variable2 of slotConnectorDataTimestepping,
  // after the already present additional field variables of the timestepping scheme.

  // The first "states" entry of statesToTransfer is the solution variable, component 0 (which is default) of the timestepping scheme and therefore
  // the timestepping scheme has already added it to the slotConnectorDataTimestepping object.
  // Now remove it because we set all connections of the CellmlAdapter here.
  slotConnectorDataTimestepping_->variable1.erase(slotConnectorDataTimestepping_->variable1.begin());

  // count number of slots that are defined in the CellmlAdapter settings
  int nSlotsCellmlAdapter = statesForTransfer_.size() + algebraicsForTransfer_.size() + parametersForTransfer_.size();


  LOG(DEBUG) << "got the following states for transfer: " << statesForTransfer_ << " (states: " << this->states()
    << "), algebraics: " << algebraicsForTransfer_ << ", parameters: " << parametersForTransfer_;

  int slotNo = 0;
  std::vector<std::string> ownSlotNames(slotNames_);
  ownSlotNames.resize(nSlotsCellmlAdapter);

  LOG(DEBUG) << "slot names of CellmlAdapter: " << slotNames_ << " -> " << ownSlotNames;

  std::vector<std::string> additionalSlotNamesTimeSteppingScheme = slotConnectorDataTimestepping_->slotNames;
  LOG(DEBUG) << "CellmlAdapter::setSlotConnectorDataTimestepping " << slotConnectorDataTimestepping_->slotNames.size()
    << " timestepping slot names that will be cleared: " << slotConnectorDataTimestepping_->slotNames << ", " << ownSlotNames.size() << " own slot names: " << ownSlotNames;

  // remove first slot name of timestepping scheme:
  if (!additionalSlotNamesTimeSteppingScheme.empty())
  {
    additionalSlotNamesTimeSteppingScheme.erase(additionalSlotNamesTimeSteppingScheme.begin());
  }

  // clear all slot names of the timestepping scheme, they will be set anew in this method
  slotConnectorDataTimestepping_->slotNames.clear();

  // ------
  // prepare states components, set field variable to global representation
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nStates>> states = this->states();

  if (!statesForTransfer_.empty())
    states->setRepresentationGlobal();

  // loop over states that should be transferred
  for (std::vector<int>::iterator iter = statesForTransfer_.begin(); iter != statesForTransfer_.end(); iter++, slotNo++)
  {
    int componentNo = *iter;

    // The state field variables have 'nStates' components and can be reused.
    std::string name = states->componentName(componentNo);
    LOG(DEBUG) << "CellmlAdapter::setSlotConnectorDataTimestepping add FieldVariable " << *states << " (" << states->name() << ") for state " << componentNo << "," << name;

    // add this component to slotConnector of data time stepping
    slotConnectorDataTimestepping_->addFieldVariable(states, componentNo);

    // add the corresponding slot name
    slotConnectorDataTimestepping_->slotNames.push_back(ownSlotNames[slotNo]);
  }

  // after all slots of "variable1" there will be the slots of the additional field variables of the timestepping scheme and then the normal slots of "variable2"

  // add slot names for the additional slots of the timestepping scheme
  slotConnectorDataTimestepping_->slotNames.insert(slotConnectorDataTimestepping_->slotNames.end(),
                                                   additionalSlotNamesTimeSteppingScheme.begin(), additionalSlotNamesTimeSteppingScheme.end());

  // ------
  // prepare algebraic components, set field variable to global representation
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nAlgebraics>> algebraics = this->algebraics();

  if (!algebraicsForTransfer_.empty())
    algebraics->setRepresentationGlobal();

  // loop over algebraics that should be transferred
  for (std::vector<int>::iterator iter = algebraicsForTransfer_.begin(); iter != algebraicsForTransfer_.end(); iter++, slotNo++)
  {
    int componentNo = *iter;

    // The algebraic and parameter field variables have 'nAlgebraics' components, but the field variables in the slotConnectorDataTimestepping_ object
    // have only 1 component. Therefore, we create new field variables with 1 components each that reuse the Petsc Vec's of the algebraic or parameter field variables.

    // get the parameters to create the new field variable
    std::string name = algebraics->componentName(componentNo);
    const std::vector<std::string> componentNames{"0"};
    const bool reuseData = true;

    // create the new field variable with only the one component, the component given by componentNo
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> newFieldVariable
      = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,1>>(*algebraics, name, componentNames, reuseData, componentNo);

    LOG(DEBUG) << "CellmlAdapterBase::setSlotConnectorDataTimestepping add FieldVariable2 (an algebraic) " << newFieldVariable << " with name " << name << " for component no " << componentNo
      << ", this reuses the data from \"" << algebraics->name() << "\".";

    // add this component to slotConnector of data time stepping
    slotConnectorDataTimestepping_->addFieldVariable2(newFieldVariable);

    // add the corresponding slot name
    slotConnectorDataTimestepping_->slotNames.push_back(ownSlotNames[slotNo]);
  }

  // ------
  // prepare parameters components, set field variable to global representation
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nAlgebraics>> parameters = this->parameters();

  if (!parametersForTransfer_.empty())
    parameters->setRepresentationGlobal();


  // add parameters components
  for (std::vector<int>::iterator iter = parametersForTransfer_.begin(); iter != parametersForTransfer_.end(); iter++, slotNo++)
  {
    int componentNo = *iter;

    // The algebraic and parameter field variables have 'nAlgebraics' components, but the field variables in the slotConnectorDataTimestepping_ object
    // have only 1 component. Therefore, we create new field variables with 1 components each that reuse the Petsc Vec's of the algebraic or paramater field variables.

    // get the parameters to create the new field variable
    std::string name = parameters->componentName(componentNo);
    const std::vector<std::string> componentNames{"0"};
    const bool reuseData = true;

    // create the new field variable with only the one component, the component given by componentNo
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> newFieldVariable
      = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,1>>(*parameters, name, componentNames, reuseData, componentNo);

    LOG(DEBUG) << "CellmlAdapterBase::setSlotConnectorDataTimestepping add FieldVariable2 (a parameter) " << newFieldVariable << " with name " << name << " for component no " << componentNo
      << ", this reuses the data from \"" << parameters->name() << "\".";

    // add this component to slotConnector of data time stepping
    slotConnectorDataTimestepping_->addFieldVariable2(newFieldVariable);

    // add the corresponding slot name
    slotConnectorDataTimestepping_->slotNames.push_back(ownSlotNames[slotNo]);
  }

  // output the slot names
  LOG(DEBUG) << "CellmlAdapterBase::setSlotConnectorData new slot names: " << slotConnectorDataTimestepping_->slotNames;
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


  if (parameterNames.size() > nAlgebraics)
  {
    LOG(FATAL) << "There can only be as many parameters as there are algebraics. This is an arbitrary restriction, if you need more parameters, try increasing the number of algebraics in the C++ source file."
      << "Now you have (" << parameterNames.size() << " parameters and the maximum possible number is " << nAlgebraics << ").";
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

//! get the parameterValues_ pointer from the parameters field variable, then the field variable can no longer be used until restoreParameterValues() gets called
template <int nStates, int nAlgebraics, typename FunctionSpaceType>
void CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
prepareParameterValues()
{
  this->parameters_->setRepresentationContiguous();
  PetscErrorCode ierr;
  Vec contiguousVec = this->parameters_->getValuesContiguous();
  ierr = VecGetArray(contiguousVec, &parameterValues_); CHKERRV(ierr);

  updateParametersPointer();
  
  // debugging output
#if 0
  PetscInt nValues;
  ierr = VecGetLocalSize(contiguousVec, &nValues); CHKERRV(ierr);

  std::stringstream s;
  for (int i = 0; i < std::min(100,nValues); i++)
    s << ",  " << i << ": " << parameterValues_[i];
  LOG(INFO) << "parameter values (field variable " << this->parameters_ << ") has " << nValues << " entries (only showing the first 100): " << s.str();
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

//! retrieve the current parameters_ field variable from the slot connector, because is could have been changed by the slot connector transfer
template <int nStates, int nAlgebraics, typename FunctionSpaceType>
void CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
updateParametersPointer()
{
  if (!slotConnectorDataTimestepping_)
  {
    LOG(DEBUG) << "updateParametersPointer: slotConnectorDataTimestepping_ is null.";
    return;
  }
  
  // loop over parameters in variable2
  int offset = algebraicsForTransfer_.size();
  
  int nInstances = this->parameters_->functionSpace()->meshPartition()->nDofsLocalWithoutGhosts();
  int nEntriesLocal = nInstances * this->parameters_->nComponents();

  LOG(DEBUG) << "updateParametersPointer: algebraicsForTransfer_.size(): " << algebraicsForTransfer_.size() 
    << ",  parametersForTransfer_.size(): " <<  parametersForTransfer_.size() << ", this->parameters_ : " << this->parameters_ 
    << ", \"" << this->parameters_->name() << "\", " << this->parameters_->partitionedPetscVec()->getCurrentRepresentationString()
    << ", nDofsLocalWithoutGhosts: " << nInstances
    << ", nComponents: " << this->parameters_->nComponents()
    << ", nEntriesLocal: " << nEntriesLocal;
  
  for (int i = 0; i < parametersForTransfer_.size(); i++)
  {
    if (offset + i < slotConnectorDataTimestepping_->variable2.size())
    {
      LOG(DEBUG) << "  slotConnectorDataTimestepping_->variable2[" << offset + i << "].values: " << slotConnectorDataTimestepping_->variable2[offset + i].values
        << ", \"" << slotConnectorDataTimestepping_->variable2[offset + i].values->name() << "\", " 
        << slotConnectorDataTimestepping_->variable2[offset + i].values->partitionedPetscVec()->getCurrentRepresentationString();
      
      // if the pointers in slotConnectorData and parameters are both set, but to different values
      if (slotConnectorDataTimestepping_->variable2[offset + i].values && this->parameters_ 
          && dynamic_cast<void*>(slotConnectorDataTimestepping_->variable2[offset + i].values.get()) != dynamic_cast<void*>(this->parameters_.get()))
      {
        LOG(DEBUG) << "CellmlAdapter::updateParametersPointer, slotConnectorDataTimestepping_->variable2 contains " << slotConnectorDataTimestepping_->variable2.size() << " entries: " 
          << offset << " algebraics for transfer and " << parametersForTransfer_.size() << " parameters for transfer."
          << " In parameter entry no. " << i << ", the field variable is " << slotConnectorDataTimestepping_->variable2[offset + i].values 
          << " (name " << slotConnectorDataTimestepping_->variable2[offset + i].values->name()  << ")"
          << " but the own stored parameters_ is " << this->parameters_ << " (name " << this->parameters_->name() 
          << "). Now set the parameters_ pointer to the field variable in slotConnectorData.";
          
        // copy values to the respective location within the contiguous parameterValues_ vector
        int parameterIndex = parametersForTransfer_[i];
        std::vector<double> values;
        slotConnectorDataTimestepping_->variable2[offset + i].values->getValuesWithoutGhosts(values);
        
        LOG(DEBUG) << "copy " << values.size() << " entries (nInstances=" << nInstances << ", nAlgebraics(=nParameters)=" << nAlgebraics << ", "
          << "nDofs: " << this->parameters_->nDofsLocalWithoutGhosts() << ","
          << this->parameters_->functionSpace()->nDofsLocalWithGhosts() << "," << this->parameters_->nDofsGlobal()
          << ", nEntriesLocal: " << nEntriesLocal << ") from field variable "
          << slotConnectorDataTimestepping_->variable2[offset + i].values->name()
          << " to the location of the parameterValues pointer (parameterIndex " << parameterIndex << "), from "
          << (void *)values.data() << " to " << (void *)(parameterValues_ + parameterIndex*nAlgebraics*sizeof(double));
        
        if (parameterIndex*nInstances + values.size()-1 > nEntriesLocal)
          LOG(ERROR) << "invalid copy operation, parameterValues is accessed at " << parameterIndex << "*" << nAlgebraics << "+" << values.size()-1
            << " (=" << parameterIndex*nInstances + (values.size()-1) << "), but nEntriesLocal: " << nEntriesLocal;

        // copy raw data from values to parameterValues_
        for (int valueIndex = 0; valueIndex < values.size(); valueIndex++)
        {
          parameterValues_[parameterIndex*nInstances + valueIndex] = values[valueIndex];
        }
        //memcpy((void *)(parameterValues_ + parameterIndex*nAlgebraics*sizeof(double)), (void *)values.data(), values.size()*sizeof(double));
      }
    }
  }
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
typename CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::FieldVariablesForOutputWriter CellmlAdapter<nStates,nAlgebraics,FunctionSpaceType>::
getFieldVariablesForOutputWriter()
{
  // these field variables will be written to output files
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> geometryField
    = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,3>>(this->functionSpace_->geometryField());

  return std::make_tuple(geometryField, algebraics_, states_, parameters_);
}

} // namespace
