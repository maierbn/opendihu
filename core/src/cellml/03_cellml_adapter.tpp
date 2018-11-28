#include "cellml/03_cellml_adapter.h"

#include <Python.h>  // has to be the first included header

#include <list>
#include <sstream>

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "utility/string_utility.h"
#include "mesh/structured_regular_fixed.h"
#include "mesh/mesh_manager.h"

//#include <libcellml>    // libcellml not used here

template<int nStates_, typename FunctionSpaceType>
CellmlAdapter<nStates_,FunctionSpaceType>::
CellmlAdapter(DihuContext context) :
  CallbackHandler<nStates_,FunctionSpaceType>(context),
  Splittable()
{
  LOG(TRACE) << "CellmlAdapter constructor";
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

  CellmlAdapterBase<nStates_,FunctionSpaceType>::initialize();
  
  // load rhs routine
  this->initializeRhsRoutine();

  this->initializeCallbackFunctions();
  
  this->outputStateIndex_ = this->specificSettings_.getOptionInt("outputStateIndex", 0, PythonUtility::NonNegative);
  this->prefactor_ = this->specificSettings_.getOptionDouble("prefactor", 1.0);
  this->internalTimeStepNo_ = 0;
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
  this->internalTimeStepNo_++;

  //PetscUtility::getVectorEntries(input, states_);
  double *states, *rates;
  PetscErrorCode ierr;
  ierr = VecGetArray(input, &states); CHKERRV(ierr);   // get r/w pointer to contiguous array of the data, VecRestoreArray() needs to be called afterwards
  ierr = VecGetArray(output, &rates); CHKERRV(ierr);

  int nStatesInput, nRates;
  ierr = VecGetSize(input, &nStatesInput); CHKERRV(ierr);
  ierr = VecGetSize(output, &nRates); CHKERRV(ierr);
  VLOG(1) << "evaluateTimesteppingRightHandSideExplicit, input nStates_: " << nStatesInput << ", output nRates: " << nRates;
  if (nStatesInput != nStates_*this->nInstances_)
  {
    LOG(ERROR) << "nStatesInput=" << nStatesInput << ", nStates_=" << nStates_ << ", nInstances=" << this->nInstances_;
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

  // get new values for parameters, call callback function of python config
  if (this->setSpecificStates_ && this->internalTimeStepNo_ % this->setSpecificStatesCallInterval_ == 0)
  {
    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;

    VLOG(1) << "call setSpecificStates, this->internalTimeStepNo_ = " << this->internalTimeStepNo_ << ", this->setSpecificStatesCallInterval_: " << this->setSpecificStatesCallInterval_;
    this->setSpecificStates_((void *)this, this->nInstances_, this->internalTimeStepNo_, currentTime, states);
  }

  //              this          STATES, RATES, WANTED,                KNOWN
  if (this->rhsRoutine_)
  {
    VLOG(1) << "call rhsRoutine_ with " << this->intermediates_.size() << " intermediates, " << this->parameters_.size() << " parameters";
    VLOG(2) << "intermediates: " << this->intermediates_ << ", parameters: " << this->parameters_;

    // call actual rhs routine from cellml code
    this->rhsRoutine_((void *)this, currentTime, states, rates, this->intermediates_.data(), this->parameters_.data());
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
}

/*
template<int nStates_, typename FunctionSpaceType>
void CellmlAdapter<nStates_,FunctionSpaceType>::
evaluateTimesteppingRightHandSideImplicit(Vec& input, Vec& output, int timeStepNo, double currentTime)
{
}
*/

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
