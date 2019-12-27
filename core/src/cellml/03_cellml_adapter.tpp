#include "cellml/03_cellml_adapter.h"

#include <Python.h>  // has to be the first included header

#include <list>
#include <sstream>

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "utility/string_utility.h"
#include "mesh/structured_regular_fixed.h"
#include "mesh/mesh_manager/mesh_manager.h"
#include "function_space/function_space.h"
#include "control/stimulation_logging.h"

//#include <libcellml>    // libcellml not used here

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
CellmlAdapter<nStates_,nIntermediates_,FunctionSpaceType>::
CellmlAdapter(DihuContext context) :
  CallbackHandler<nStates_,nIntermediates_,FunctionSpaceType>(context),
  Splittable()
{
  LOG(TRACE) << "CellmlAdapter constructor";

  // initialize filename in stimulation logging class from current settings
  Control::StimulationLogging logging(this->specificSettings_);
}

//! constructor from other CellmlAdapter, preserves the outputManager and context
template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
CellmlAdapter<nStates_,nIntermediates_,FunctionSpaceType>::
CellmlAdapter(const CellmlAdapter &rhs, std::shared_ptr<FunctionSpace> functionSpace) :
  CallbackHandler<nStates_,nIntermediates_,FunctionSpaceType>(rhs.context_, true),
  Splittable()
{
  LOG(TRACE) << "CellmlAdapter constructor from rhs";

  this->functionSpace_ = functionSpace;
  this->outputWriterManager_ = rhs.outputWriterManager_;

  return;

  // copy member variables from rhs
  this->specificSettings_ = rhs.specificSettings_;
  this->setParameters_ = rhs.setParameters_;
  this->setSpecificParameters_ = rhs.setSpecificParameters_;
  this->setSpecificStates_ = rhs.setSpecificStates_;
  this->handleResult_ = rhs.handleResult_;
  this->pythonSetParametersFunction_ = rhs.pythonSetParametersFunction_;
  this->pythonSetSpecificParametersFunction_ = rhs.pythonSetSpecificParametersFunction_;
  this->pythonSetSpecificStatesFunction_ = rhs.pythonSetSpecificStatesFunction_;
  this->pythonHandleResultFunction_ = rhs.pythonHandleResultFunction_;
  this->pySetFunctionAdditionalParameter_ = rhs.pySetFunctionAdditionalParameter_;
  this->pyHandleResultFunctionAdditionalParameter_ = rhs.pyHandleResultFunctionAdditionalParameter_;
  this->pyGlobalNaturalDofsList_ = rhs.pyGlobalNaturalDofsList_;

  this->nInstances_ = this->functionSpace_->nNodesLocalWithoutGhosts();
  assert(this->nInstances_ > 1);

  // copy member variables from rhs
  this->parametersUsedAsIntermediate_ = rhs.parametersUsedAsIntermediate_;
  this->parametersUsedAsConstant_ = rhs.parametersUsedAsConstant_;
  this->stateNames_ = rhs.stateNames_;

  this->sourceFilename_ = rhs.sourceFilename_;
  this->nIntermediates_ = rhs.nIntermediates_;
  this->nIntermediatesInSource_ = rhs.nIntermediatesInSource_;
  this->nParameters_ = rhs.nParameters_;
  this->nConstants_ = rhs.nConstants_;
  this->internalTimeStepNo_ = rhs.internalTimeStepNo_;
  this->inputFileTypeOpenCMISS_ = rhs.inputFileTypeOpenCMISS_;

  // allocate data vectors
  //this->intermediates_.resize(this->nIntermediates_*this->nInstances_);
  this->parameters_.resize(this->nParameters_*this->nInstances_);

  // copy rhs parameter values to parameters, it is assumed that the parameters are the same for every instance
  for (int instanceNo = 0; instanceNo < this->nInstances_; instanceNo++)
  {
    for (int j = 0; j < this->nParameters_; j++)
    {
      this->parameters_[j*this->nInstances_ + instanceNo] = rhs.parameters_[j];
    }
  }

  LOG(DEBUG) << "Initialize CellML with nInstances = " << this->nInstances_ << ", nParameters_ = " << this->nParameters_
    << ", nStates = " << nStates << ", nIntermediates = " << this->nIntermediates();
}

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
constexpr int CellmlAdapter<nStates_,nIntermediates_,FunctionSpaceType>::
nStates()
{
  return nStates_;
}

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
void CellmlAdapter<nStates_,nIntermediates_,FunctionSpaceType>::
reset()
{
  this->internalTimeStepNo_ = 0;
}
  
template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
void CellmlAdapter<nStates_,nIntermediates_,FunctionSpaceType>::
initialize()
{
  LOG(TRACE) << "CellmlAdapter<nStates_,nIntermediates_,FunctionSpaceType>::initialize";

  Control::PerformanceMeasurement::start("durationInitCellml");

  CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::initialize();
  
  // load rhs routine
  this->initializeRhsRoutine();

  this->initializeCallbackFunctions();
  
  Control::PerformanceMeasurement::stop("durationInitCellml");

  if (this->specificSettings_.hasKey("prefactor"))
  {
    LOG(WARNING) << this->specificSettings_ << "[\"prefactor\"] is no longer an option! There is no more functionality to scale values during transfer.";
  }

  this->internalTimeStepNo_ = 0;

  this->setSpecificStatesCallFrequency_ = this->specificSettings_.getOptionDouble("setSpecificStatesCallFrequency", 0.0);
  this->specificSettings_.getOptionVector("setSpecificStatesFrequencyJitter", this->setSpecificStatesFrequencyJitter_);
  this->setSpecificStatesRepeatAfterFirstCall_ = this->specificSettings_.getOptionDouble("setSpecificStatesRepeatAfterFirstCall", 0.0);
  this->setSpecificStatesCallEnableBegin_ = this->specificSettings_.getOptionDouble("setSpecificStatesCallEnableBegin", 0.0);

  // initialize the lastCallSpecificStatesTime_
  this->currentJitter_ = 0;
  this->jitterIndex_ = 0;
  this->lastCallSpecificStatesTime_ = this->setSpecificStatesCallEnableBegin_ - 1e-13 - 1./(this->setSpecificStatesCallFrequency_+this->currentJitter_);
}

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
void CellmlAdapter<nStates_,nIntermediates_,FunctionSpaceType>::
initializeForImplicitTimeStepping()
{
}

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
void CellmlAdapter<nStates_,nIntermediates_,FunctionSpaceType>::
setRankSubset(Partition::RankSubset rankSubset)
{
  // do nothing because we don't have stored data here (the data on which the computation is performed comes in evaluateTimesteppingRightHandSide from parameters) 
}

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
void CellmlAdapter<nStates_,nIntermediates_,FunctionSpaceType>::
evaluateTimesteppingRightHandSideExplicit(Vec& input, Vec& output, int timeStepNo, double currentTime)
{
  // get raw pointers from Petsc data structures
  //PetscUtility::getVectorEntries(input, states_);
  double *states, *rates;
  double *intermediatesData;
  PetscErrorCode ierr;
  ierr = VecGetArray(input, &states); CHKERRV(ierr);   // get r/w pointer to contiguous array of the data, VecRestoreArray() needs to be called afterwards
  ierr = VecGetArray(output, &rates); CHKERRV(ierr);
  ierr = VecGetArray(this->data_.intermediates()->getValuesContiguous(), &intermediatesData); CHKERRV(ierr);

  // get sizes of input and output Vecs
  int nStatesInput, nRates, nIntermediates = 101;
  ierr = VecGetSize(input, &nStatesInput); CHKERRV(ierr);
  ierr = VecGetSize(output, &nRates); CHKERRV(ierr);
  ierr = VecGetLocalSize(this->data_.intermediates()->getValuesContiguous(), &nIntermediates); CHKERRV(ierr);

  //double intermediatesData[101];

  VLOG(1) << "intermediates array has " << nIntermediates << " entries";
  nIntermediates = nIntermediates/this->nInstances_;

  // check validity of sizes
  VLOG(1) << "evaluateTimesteppingRightHandSideExplicit, input nStates_: " << nStatesInput << ", output nRates: " << nRates;
  VLOG(1) << "timeStepNo: " << timeStepNo << ", currentTime: " << currentTime << ", internalTimeStepNo: " << this->internalTimeStepNo_;

  if (nStatesInput != nStates_*this->nInstances_)
  {
    LOG(ERROR) << "nStatesInput does not match nStates and nInstances! nStatesInput=" << nStatesInput << ", nStates_=" << nStates_ << ", nInstances=" << this->nInstances_;
  }
  assert (nStatesInput == nStates_*this->nInstances_);
  assert (nRates == nStates_*this->nInstances_);
  if (nIntermediates != nIntermediates_)
  {
    LOG(FATAL) << "nInstances: " << this->nInstances_ << ", nIntermediates (size of vector / nInstances): " << nIntermediates << ", nIntermediates_: " << nIntermediates_;
  }
  assert (nIntermediates == nIntermediates_);

  //LOG(DEBUG) << " evaluateTimesteppingRightHandSide: nInstances=" << this->nInstances_ << ", nStates_=" << nStates_;
  
  // get new values for parameters, call callback function of python config
  if (this->setParameters_ && this->internalTimeStepNo_ % this->setParametersCallInterval_ == 0)
  {
    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;
    
    VLOG(1) << "call setParameters";
    this->setParameters_((void *)this, this->nInstances_, this->internalTimeStepNo_, currentTime, this->parameters_);
  }

  // get new values for parameters, call callback function of python config
  if (this->setSpecificParameters_ && this->internalTimeStepNo_ % this->setSpecificParametersCallInterval_ == 0)
  {
    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;

    VLOG(1) << "call setSpecificParameters";
    this->setSpecificParameters_((void *)this, this->nInstances_, this->internalTimeStepNo_, currentTime, this->parameters_);
  }

  VLOG(1) << "currentTime: " << currentTime << ", lastCallSpecificStatesTime_: " << this->lastCallSpecificStatesTime_
    << ", setSpecificStatesCallFrequency_: " << this->setSpecificStatesCallFrequency_ << ", "
    << this->lastCallSpecificStatesTime_ + 1./this->setSpecificStatesCallFrequency_;
  VLOG(1) << "this->setSpecificStates_? " << (this->setSpecificStates_? "true" : "false")
    << ", this->setSpecificStatesCallInterval_: " << this->setSpecificStatesCallInterval_ << " != 0? " << (this->setSpecificStatesCallInterval_ != 0? "true" : "false")
    << ", (this->internalTimeStepNo_ % this->setSpecificStatesCallInterval_) = " << this->internalTimeStepNo_  << " % " << this->setSpecificStatesCallInterval_
    << ", this->setSpecificStatesCallFrequency_= " << this->setSpecificStatesCallFrequency_ << " != 0.0? " << (this->setSpecificStatesCallFrequency_ != 0.0? "true" : "false")
    << ", currentTime=" << currentTime << " >= " << this->lastCallSpecificStatesTime_ << " + " << 1./this->setSpecificStatesCallFrequency_ << " = " << this->lastCallSpecificStatesTime_ + 1./this->setSpecificStatesCallFrequency_ << "? "
    << (currentTime >= this->lastCallSpecificStatesTime_ + 1./this->setSpecificStatesCallFrequency_? "true" : "false");

  bool stimulate = false;

  // get new values for parameters, call callback function of python config
  if (this->setSpecificStates_
      && (
          (this->setSpecificStatesCallInterval_ != 0 && this->internalTimeStepNo_ % this->setSpecificStatesCallInterval_ == 0)
          || (this->setSpecificStatesCallFrequency_ != 0.0 && currentTime >= this->lastCallSpecificStatesTime_ + 1./(this->setSpecificStatesCallFrequency_+this->currentJitter_)
              && currentTime >= this->setSpecificStatesCallEnableBegin_-1e-13)
         )
     )
  {
    stimulate = true;

    // if current stimulation is over
    if (currentTime - (this->lastCallSpecificStatesTime_ + 1./(this->setSpecificStatesCallFrequency_+this->currentJitter_)) > this->setSpecificStatesRepeatAfterFirstCall_)
    {
      // advance time of last call to specificStates
      this->lastCallSpecificStatesTime_ += 1./(this->setSpecificStatesCallFrequency_+this->currentJitter_);

      // get new jitter value
      double jitterFactor = this->setSpecificStatesFrequencyJitter_[this->jitterIndex_ % this->setSpecificStatesFrequencyJitter_.size()];
      this->currentJitter_ = jitterFactor * this->setSpecificStatesCallFrequency_;
      this->jitterIndex_++;

      stimulate = false;
    }
    VLOG(1) << "next call to setSpecificStates,this->setSpecificStatesCallFrequency_: " << this->setSpecificStatesCallFrequency_ << ", set lastCallSpecificStatesTime_ to " << this->lastCallSpecificStatesTime_;

    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;
  }

  static bool currentlyStimulating = false;
  if (stimulate)
  {
    // if this is the first point in time of the current stimulation, log stimulation time
    if (!currentlyStimulating)
    {
      currentlyStimulating = true;
      int fiberNoGlobal = -1;
      if (this->pySetFunctionAdditionalParameter_)
      {
        fiberNoGlobal = PythonUtility::convertFromPython<int>::get(this->pySetFunctionAdditionalParameter_);
      }
      Control::StimulationLogging::logStimulationBegin(currentTime, -1, fiberNoGlobal);
    }

    VLOG(1) << "call setSpecificStates, this->internalTimeStepNo_ = " << this->internalTimeStepNo_ << ", this->setSpecificStatesCallInterval_: " << this->setSpecificStatesCallInterval_;
    VLOG(1) << "currentTime: " << currentTime << ", call setSpecificStates, this->internalTimeStepNo_ = " << this->internalTimeStepNo_ << ", this->setSpecificStatesCallInterval_: " << this->setSpecificStatesCallInterval_;
    this->setSpecificStates_((void *)this, this->nInstances_, this->internalTimeStepNo_, currentTime, states);
  }
  else
  {
    currentlyStimulating = false;
  }

  // call actual rhs method
  //              this          STATES, RATES, WANTED,                KNOWN
  if (this->rhsRoutine_)
  {
    VLOG(1) << "call rhsRoutine_ with " << nIntermediates << " intermediates, " << this->parameters_.size() << " parameters";
    VLOG(2) << "parameters: " << this->parameters_;

    //Control::PerformanceMeasurement::start("rhsEvaluationTime");  // commented out because it takes too long in this very inner loop
    // call actual rhs routine from cellml code
    this->rhsRoutine_((void *)this, currentTime, states, rates, intermediatesData, this->parameters_.data());
    //Control::PerformanceMeasurement::stop("rhsEvaluationTime");
  }

  // handle intermediates, call callback function of python config
  if (this->handleResult_ && this->internalTimeStepNo_ % this->handleResultCallInterval_ == 0)
  {
    int nStatesInput;
    VecGetSize(input, &nStatesInput);

    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;
    
    VLOG(1) << "call handleResult with in total " << nStatesInput << " states, " << nIntermediates << " intermediates";
    this->handleResult_((void *)this, this->nInstances_, this->internalTimeStepNo_, currentTime, states, intermediatesData);
  }

  //PetscUtility::setVector(rates_, output);
  // give control of data back to Petsc
  ierr = VecRestoreArray(input, &states); CHKERRV(ierr);
  ierr = VecRestoreArray(output, &rates); CHKERRV(ierr);
  ierr = VecRestoreArray(this->data_.intermediates()->getValuesContiguous(), &intermediatesData); CHKERRV(ierr);

  VLOG(1) << "intermediates: " << *this->data_.intermediates();

  this->internalTimeStepNo_++;
}

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
void CellmlAdapter<nStates_,nIntermediates_,FunctionSpaceType>::
getComponentNames(std::vector<std::string> &stateNames)
{
  this->getStateNames(stateNames);
}

//! return the mesh
template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
std::shared_ptr<FunctionSpaceType> CellmlAdapter<nStates_,nIntermediates_,FunctionSpaceType>::
functionSpace()
{
  return CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::functionSpace();
}

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
template<typename FunctionSpaceType2>
bool CellmlAdapter<nStates_,nIntermediates_,FunctionSpaceType>::
setInitialValues(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nStates_>> initialValues)
{
  return CellmlAdapterBase<nStates_,nIntermediates_,FunctionSpaceType>::template setInitialValues<FunctionSpaceType2>(initialValues);
}
