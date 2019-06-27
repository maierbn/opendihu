#include "specialized_solver/reaction_diffusion_accelerator.h"

#include <Python.h>  // has to be the first included header

namespace TimeSteppingScheme
{

template<int nStates, int nIntermediates>
ReactionDiffusionAccelerator<nStates,nIntermediates>::
ReactionDiffusionAccelerator(DihuContext context) :
  context_(context["ReactionDiffusionAccelerator"]),
  data_(this->context_), initialized_(false)
{
  // get python config
  this->specificSettings_ = this->context_.getPythonConfig();

  // read in the durationLogKey
  if (specificSettings_.hasKey("durationLogKey"))
  {
    this->durationLogKey_ = specificSettings_.getOptionString("durationLogKey", "");
  }

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_, this->specificSettings_);
}

template<int nStates, int nIntermediates>
bool ReactionDiffusionAccelerator<nStates,nIntermediates>::
setInitialValues(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> initialValues)
{
  CellMLAdapter<nStates,nIntermediates,FunctionSpaceType>::setInitialValues(initialValues);
  return true;
}

template<int nStates, int nIntermediates>
void ReactionDiffusionAccelerator<nStates,nIntermediates>::
advanceTimeSpan()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);

  // write current output values
  this->outputWriterManager_.writeOutput(this->data_, 0, endTime_);
}

template<int nStates, int nIntermediates>
void ReactionDiffusionAccelerator<nStates,nIntermediates>::
run()
{
  // initialize everything
  initialize();

  this->advanceTimeSpan();
}

template<int nStates, int nIntermediates>
void ReactionDiffusionAccelerator<nStates,nIntermediates>::
initialize()
{
  if (this->initialized_)
    return;
  TimeSteppingScheme::initialize();
  
  assert(this->specificSettings_.pyObject());
  
  std::shared_ptr<FunctionSpaceType> functionSpace = this->context_.meshManager()->functionSpace<FunctionSpaceType>(this->specificSettings_);
 
  // initialize the data object
  data_.setFunctionSpace(functionSpace);
  data_.initialize();

  LOG(DEBUG) << "initialization done";
  this->initialized_ = true;
}

template<int nStates, int nIntermediates>
void ReactionDiffusionAccelerator<nStates,nIntermediates>::reset()
{
  this->initialized_ = false;
}

//! return whether the underlying discretizableInTime object has a specified mesh type and is not independent of the mesh type
template<int nStates, int nIntermediates>
bool ReactionDiffusionAccelerator<nStates,nIntermediates>::
knowsMeshType()
{
  return true;
}

template<int nStates, int nIntermediates>
typename ReactionDiffusionAccelerator<nStates,nIntermediates>::Data &ReactionDiffusionAccelerator<nStates,nIntermediates>::
data()
{
  return data_;
}

//! get the data that will be transferred in the operator splitting to the other term of the splitting
//! the transfer is done by the solution_vector_mapping class
template<int nStates, int nIntermediates>
typename ReactionDiffusionAccelerator<nStates,nIntermediates>::TransferableSolutionDataType
ReactionDiffusionAccelerator<nStates,nIntermediates>::
getSolutionForTransfer()
{
  return this->data_.intermediates();
}

//! output the given data for debugging
template<int nStates, int nIntermediates>
std::string ReactionDiffusionAccelerator<nStates,nIntermediates>::
getString(typename ReactionDiffusionAccelerator<nStates,nIntermediates>::TransferableSolutionDataType &data)
{
  std::stringstream s;
  s << "<reaction_diffusion_accelerator:" << *data << ">";
  return s.str();
}

} // namespace TimeSteppingScheme
