#include "cellml/03_cellml_adapter.h"

#include <Python.h>  // has to be the first included header

#include <list>
#include <sstream>

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "utility/string_utility.h"
#include "mesh/structured_regular_fixed.h"
#include "mesh/mesh_manager.h"
#include "function_space/function_space.h"

//#include <libcellml>    // libcellml not used here

template<int nStates_, typename FunctionSpaceType>
CellmlAdapter<nStates_,FunctionSpaceType>::
CellmlAdapter(DihuContext context) :
  CallbackHandler<nStates_,FunctionSpaceType>(context),
  Splittable()
{
  LOG(TRACE) << "CellmlAdapter constructor";
}

//! constructor from other CellmlAdapter, preserves the outputManager and context
template<int nStates_, typename FunctionSpaceType>
CellmlAdapter<nStates_,FunctionSpaceType>::
CellmlAdapter(const CellmlAdapter &rhs, std::shared_ptr<FunctionSpace> functionSpace) :
  CallbackHandler<nStates_,FunctionSpaceType>(rhs.context_, true),
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
  this->nParameters_ = rhs.nParameters_;
  this->nConstants_ = rhs.nConstants_;
  this->outputStateIndex_ = rhs.outputStateIndex_;
  this->prefactor_ = rhs.prefactor_;
  this->internalTimeStepNo_ = rhs.internalTimeStepNo_;
  this->inputFileTypeOpenCMISS_ = rhs.inputFileTypeOpenCMISS_;

  // allocate data vectors
  this->intermediates_.resize(this->nIntermediates_*this->nInstances_);
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
    << ", nStates = " << nStates << ", nIntermediates = " << this->nIntermediates_;
}

template<int nStates_, typename FunctionSpaceType>
constexpr int CellmlAdapter<nStates_,FunctionSpaceType>::
nStates()
{
  return nStates_;
}

template<int nStates_, typename FunctionSpaceType>
void CellmlAdapter<nStates_,FunctionSpaceType>::
reset()
{
  this->internalTimeStepNo_ = 0;
}
  
template<int nStates_, typename FunctionSpaceType>
void CellmlAdapter<nStates_,FunctionSpaceType>::
initialize()
{
  LOG(TRACE) << "CellmlAdapter<nStates_,FunctionSpaceType>::initialize";

  Control::PerformanceMeasurement::start("durationInitCellml");

  CellmlAdapterBase<nStates_,FunctionSpaceType>::initialize();
  
  // load rhs routine
  this->initializeRhsRoutine();

  this->initializeCallbackFunctions();
  
  Control::PerformanceMeasurement::stop("durationInitCellml");

  this->outputStateIndex_ = this->specificSettings_.getOptionInt("outputStateIndex", 0, PythonUtility::NonNegative);
  this->prefactor_ = this->specificSettings_.getOptionDouble("prefactor", 1.0);

  this->internalTimeStepNo_ = 0;

  this->setSpecificStatesCallFrequency_ = this->specificSettings_.getOptionDouble("setSpecificStatesCallFrequency", 0.0);

  // initialize the lastCallSpecificStatesTime_ to something negative, such that the condition is fullfilled already in the fisrrt iteration
  this->lastCallSpecificStatesTime_ = -2*this->setSpecificStatesCallFrequency_;

}

template<int nStates_, typename FunctionSpaceType>
void CellmlAdapter<nStates_,FunctionSpaceType>::
initializeForImplicitTimeStepping()
{
}

template<int nStates_, typename FunctionSpaceType>
void CellmlAdapter<nStates_,FunctionSpaceType>::
setRankSubset(Partition::RankSubset rankSubset)
{
  // do nothing because we don't have stored data here (the data on which the computation is performed comes in evaluateTimesteppingRightHandSide from parameters) 
}

template<int nStates_, typename FunctionSpaceType>
void CellmlAdapter<nStates_,FunctionSpaceType>::
evaluateTimesteppingRightHandSideExplicit(Vec& input, Vec& output, int timeStepNo, double currentTime)
{
  //PetscUtility::getVectorEntries(input, states_);
  double *states, *rates;
  PetscErrorCode ierr;
  ierr = VecGetArray(input, &states); CHKERRV(ierr);   // get r/w pointer to contiguous array of the data, VecRestoreArray() needs to be called afterwards
  ierr = VecGetArray(output, &rates); CHKERRV(ierr);

  int nStatesInput, nRates;
  ierr = VecGetSize(input, &nStatesInput); CHKERRV(ierr);
  ierr = VecGetSize(output, &nRates); CHKERRV(ierr);

  VLOG(1) << "evaluateTimesteppingRightHandSideExplicit, input nStates_: " << nStatesInput << ", output nRates: " << nRates;
  VLOG(1) << "timeStepNo: " << timeStepNo << ", currentTime: " << currentTime << ", internalTimeStepNo: " << this->internalTimeStepNo_;

  if (nStatesInput != nStates_*this->nInstances_)
  {
    LOG(ERROR) << "nStatesInput does not match nStates and nInstances! nStatesInput=" << nStatesInput << ", nStates_=" << nStates_ << ", nInstances=" << this->nInstances_;
  }
  assert (nStatesInput == nStates_*this->nInstances_);
  assert (nRates == nStates_*this->nInstances_);

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

  LOG(DEBUG) << "currentTime: " << currentTime << ", lastCallSpecificStatesTime_: " << this->lastCallSpecificStatesTime_
    << ", setSpecificStatesCallFrequency_: " << this->setSpecificStatesCallFrequency_ << ", "
    << this->lastCallSpecificStatesTime_ + 1./this->setSpecificStatesCallFrequency_;

  // get new values for parameters, call callback function of python config
  if (this->setSpecificStates_
      && (
          (this->setSpecificStatesCallInterval_ != 0 && this->internalTimeStepNo_ % this->setSpecificStatesCallInterval_ == 0)
          || (this->setSpecificStatesCallFrequency_ != 0.0 && currentTime >= this->lastCallSpecificStatesTime_ + 1./this->setSpecificStatesCallFrequency_)
         )
     )
  {
    if (this->lastCallSpecificStatesTime_ < 0)
    {
      LOG(DEBUG) << "initial call to setSpecificStates, set lastCallSpecificStatesTime_ to 0, was " << this->lastCallSpecificStatesTime_;
      this->lastCallSpecificStatesTime_ = 0;
    }
    else
    {
      this->lastCallSpecificStatesTime_ += 1./this->setSpecificStatesCallFrequency_;
      LOG(DEBUG) << "next call to setSpecificStates,this->setSpecificStatesCallFrequency_: " << this->setSpecificStatesCallFrequency_ << ", set lastCallSpecificStatesTime_ to " << this->lastCallSpecificStatesTime_;
    }

    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;

    LOG(DEBUG) << "call setSpecificStates, this->internalTimeStepNo_ = " << this->internalTimeStepNo_ << ", this->setSpecificStatesCallInterval_: " << this->setSpecificStatesCallInterval_;
    LOG(DEBUG) << "currentTime: " << currentTime << ", call setSpecificStates, this->internalTimeStepNo_ = " << this->internalTimeStepNo_ << ", this->setSpecificStatesCallInterval_: " << this->setSpecificStatesCallInterval_;
    this->setSpecificStates_((void *)this, this->nInstances_, this->internalTimeStepNo_, currentTime, states);
  }

  //              this          STATES, RATES, WANTED,                KNOWN
  if (this->rhsRoutine_)
  {
    VLOG(1) << "call rhsRoutine_ with " << this->intermediates_.size() << " intermediates, " << this->parameters_.size() << " parameters";
    VLOG(2) << "intermediates: " << this->intermediates_ << ", parameters: " << this->parameters_;

    //Control::PerformanceMeasurement::start("rhsEvaluationTime");  // commented out because it takes too long in this very inner loop
    // call actual rhs routine from cellml code
    this->rhsRoutine_((void *)this, currentTime, states, rates, this->intermediates_.data(), this->parameters_.data());
    //Control::PerformanceMeasurement::stop("rhsEvaluationTime");
  }

  // handle intermediates, call callback function of python config
  if (this->handleResult_ && this->internalTimeStepNo_ % this->handleResultCallInterval_ == 0)
  {
    int nStatesInput;
    VecGetSize(input, &nStatesInput);

    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;
    
    VLOG(1) << "call handleResult with in total " << nStatesInput << " states, " << this->intermediates_.size() << " intermediates";
    this->handleResult_((void *)this, this->nInstances_, this->internalTimeStepNo_, currentTime, states, this->intermediates_.data());
  }

  //PetscUtility::setVector(rates_, output);
  // give control of data back to Petsc
  ierr = VecRestoreArray(input, &states); CHKERRV(ierr);
  ierr = VecRestoreArray(output, &rates); CHKERRV(ierr);

  this->internalTimeStepNo_++;
}

//! return false because the object is independent of mesh type
template<int nStates_, typename FunctionSpaceType>
bool CellmlAdapter<nStates_,FunctionSpaceType>::
knowsMeshType()
{
  return CellmlAdapterBase<nStates_,FunctionSpaceType>::knowsMeshType();
}

template<int nStates_, typename FunctionSpaceType>
void CellmlAdapter<nStates_,FunctionSpaceType>::
getComponentNames(std::vector<std::string> &stateNames)
{
  this->getStateNames(stateNames);
}

//! return the mesh
template<int nStates_, typename FunctionSpaceType>
std::shared_ptr<FunctionSpaceType> CellmlAdapter<nStates_,FunctionSpaceType>::
functionSpace()
{
  return CellmlAdapterBase<nStates_,FunctionSpaceType>::functionSpace();
}

template<int nStates_, typename FunctionSpaceType>
template<typename FunctionSpaceType2>
bool CellmlAdapter<nStates_,FunctionSpaceType>::
setInitialValues(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nStates_>> initialValues)
{
  return CellmlAdapterBase<nStates_,FunctionSpaceType>::template setInitialValues<FunctionSpaceType2>(initialValues);
}
