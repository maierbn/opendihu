#include "time_stepping_scheme/02_time_stepping_scheme_ode.h"

#include <Python.h>  // has to be the first included header
#include <vector>

#include "utility/python_utility.h"

namespace TimeSteppingScheme
{

template<typename FunctionSpaceType, int nComponents>
TimeSteppingSchemeOdeBase<FunctionSpaceType,nComponents>::
TimeSteppingSchemeOdeBase(DihuContext context, std::string name) :
TimeSteppingScheme(context[name]), initialized_(false), name_(name)
{
  // get python config
  this->specificSettings_ = this->context_.getPythonConfig();

  checkForNanInf_ = this->specificSettings_.getOptionBool("checkForNanInf", false);

  // initialize output writers
  this->outputWriterManager_.initialize(context_, this->specificSettings_);
}

template<typename FunctionSpaceType, int nComponents>
typename Data::TimeStepping<FunctionSpaceType, nComponents> &TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
data()
{
  return *data_;
}

template<typename FunctionSpaceType, int nComponents>
void TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
setInitialValues()
{
  // set initial values as given in settings, or set to zero if not given
  std::vector<double> localValues;

  bool inputMeshIsGlobal = this->specificSettings_.getOptionBool("inputMeshIsGlobal", true);
  if (inputMeshIsGlobal)
  {
    assert(this->data_);
    assert(this->data_->functionSpace());
    const int nDofsGlobal = this->data_->functionSpace()->nDofsGlobal();
    LOG(DEBUG) << "setInitialValues, nDofsGlobal = " << nDofsGlobal;

    this->specificSettings_.getOptionVector("initialValues", nDofsGlobal, localValues);

    this->data_->functionSpace()->meshPartition()->extractLocalDofsWithoutGhosts(localValues);
  }
  else
  {
    const int nDofsLocal = this->data_->functionSpace()->nDofsLocalWithoutGhosts();
    this->specificSettings_.getOptionVector("initialValues", nDofsLocal, localValues);
  }
  VLOG(1) << "set initial values to " << localValues;

  // set the first component of the solution variable by the given values
  data_->solution()->setValuesWithoutGhosts(0, localValues);

  VLOG(1) << data_->solution();
}

template<typename FunctionSpaceType, int nComponents>
void TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
setRankSubset(Partition::RankSubset rankSubset)
{
  data_->setRankSubset(rankSubset);
}

template<typename FunctionSpaceType, int nComponents>
void TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
reset()
{
  TimeSteppingScheme::reset();
  initialized_ = false;

  this->outputWriterManager_.initialize(context_, this->specificSettings_);
}

template<typename FunctionSpaceType, int nComponents>
void TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
initialize()
{
  if (initialized_)
    return;

  // initialize the parent class
  TimeSteppingScheme::initialize();

  initialized_ = true;
}

template<typename FunctionSpaceType, int nComponents>
void TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
run()
{
  // initialize
  this->initialize();

  // do simulations
  this->advanceTimeSpan();
}

template<typename FunctionSpaceType, int nComponents>
void TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
checkForNanInf(int timeStepNo, double currentTime)
{
  if (checkForNanInf_)
  {
    if (this->data_->solution()->containsNanOrInf())
    {
      LOG(ERROR) << "In " << name_ << ", timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime << ": Solution contains Nan or Inf. "
        << "This probably means that the timestep width, " << this->timeStepWidth_ << ", is too high.";
      LOG(ERROR) << *this->data_->solution();
      LOG(FATAL) << "Abort because of nan or inf in solution. Set option \"checkForNanInf\": False to avoid this.";
    }
  }
}

template<typename FunctionSpaceType, int nComponents>
std::shared_ptr<typename TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::OutputConnectorDataType>
TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
getOutputConnectorData()
{
  //LOG(TRACE) << "\ncall prepareForGetOutputConnectorData";
  prepareForGetOutputConnectorData();

  return this->data_->getOutputConnectorData();
}


} // namespace
