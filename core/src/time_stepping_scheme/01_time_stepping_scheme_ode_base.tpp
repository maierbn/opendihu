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
  std::vector<VecD<nComponents>> localValues;

  if (this->specificSettings_.hasKey("initialValues") && !this->specificSettings_.isEmpty("initialValues"))
  {

    // determine if the initial values are given as global or local array
    bool inputMeshIsGlobal = this->specificSettings_.getOptionBool("inputMeshIsGlobal", true);
    if (inputMeshIsGlobal)
    {
      // if the settings specify a global list of values, extract the local values
      assert(this->data_);
      assert(this->data_->functionSpace());

      // get number of global dofs, i.e. number of values in global list
      const int nDofsGlobal = this->data_->functionSpace()->nDofsGlobal();
      LOG(DEBUG) << "setInitialValues, nDofsGlobal = " << nDofsGlobal;

      this->specificSettings_.getOptionVector("initialValues", nDofsGlobal, localValues);

      // extract only the local dofs out of the list of global values
      this->data_->functionSpace()->meshPartition()->extractLocalDofsWithoutGhosts(localValues);
    }
    else
    {
      // input is already only the local dofs, use all
      const int nDofsLocal = this->data_->functionSpace()->nDofsLocalWithoutGhosts();
      this->specificSettings_.getOptionVector("initialValues", nDofsLocal, localValues);
    }
    VLOG(1) << "set initial values to " << localValues;

    // set the first component of the solution variable by the given values
    this->data_->solution()->setValuesWithoutGhosts(localValues);

    VLOG(1) << this->data_->solution();
  }
  else
  {
    this->data_->solution()->zeroEntries();
  }
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

//! call the output writer on the data object, output files will contain currentTime, with callCountIncrement !=1 output timesteps can be skipped
template<typename FunctionSpaceType, int nComponents>
void TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
callOutputWriter(int timeStepNo, double currentTime, int callCountIncrement)
{
  this->outputWriterManager_.writeOutput(*this->data_, timeStepNo, currentTime, callCountIncrement);
}

template<typename FunctionSpaceType, int nComponents>
void TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
checkForNanInf(int timeStepNo, double currentTime)
{
  if (checkForNanInf_)
  {
    // only check every 10th time step
    static int counter = 0;
    if (counter % 10 == 0)
    {
      if (this->data_->solution()->containsNanOrInf())
      {
        LOG(ERROR) << "In " << name_ << ", timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime << ": Solution contains Nan or Inf. "
          << "This probably means that the timestep width, " << this->timeStepWidth_ << ", is too high.";
        LOG(ERROR) << *this->data_->solution();
        LOG(FATAL) << "Abort because of nan or inf in solution. Set option \"checkForNanInf\": False to avoid this.";
      }
    }
    counter++;
  }
}

template<typename FunctionSpaceType, int nComponents>
std::shared_ptr<typename TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::SlotConnectorDataType>
TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::
getSlotConnectorData()
{
  prepareForGetSlotConnectorData();

  return this->data_->getSlotConnectorData();
}


} // namespace
