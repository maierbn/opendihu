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
#include "control/diagnostic_tool/stimulation_logging.h"

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
  PetscInt nStatesInput, nRates, nIntermediates = 101;
  ierr = VecGetSize(input, &nStatesInput); CHKERRV(ierr);
  ierr = VecGetSize(output, &nRates); CHKERRV(ierr);
  ierr = VecGetLocalSize(this->data_.intermediates()->getValuesContiguous(), &nIntermediates); CHKERRV(ierr);

  // get parameterValues_ vector
  this->data_.prepareParameterValues();

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
    this->setParameters_((void *)this, this->nInstances_, this->internalTimeStepNo_, currentTime, this->data_.parameterValues(), this->cellmlSourceCodeGenerator_.nParameters());
  }

  // get new values for parameters, call callback function of python config
  if (this->setSpecificParameters_ && this->internalTimeStepNo_ % this->setSpecificParametersCallInterval_ == 0)
  {
    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;

    VLOG(1) << "call setSpecificParameters";
    this->setSpecificParameters_((void *)this, this->nInstances_, this->internalTimeStepNo_, currentTime, this->data_.parameterValues(), this->cellmlSourceCodeGenerator_.nParameters());
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
      double jitterFactor = 0.0;
      if (this->setSpecificStatesFrequencyJitter_.size() > 0)
        jitterFactor = this->setSpecificStatesFrequencyJitter_[this->jitterIndex_ % this->setSpecificStatesFrequencyJitter_.size()];
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
    VLOG(1) << "call rhsRoutine_ with " << nIntermediates << " intermediates";

    //Control::PerformanceMeasurement::start("rhsEvaluationTime");  // commented out because it takes too long in this very inner loop
    // call actual rhs routine from cellml code
    this->rhsRoutine_((void *)this, currentTime, states, rates, intermediatesData, this->data_.parameterValues());
    //Control::PerformanceMeasurement::stop("rhsEvaluationTime");
  }

  // handle intermediates, call callback function of python config
  if (this->handleResult_ && this->internalTimeStepNo_ % this->handleResultCallInterval_ == 0)
  {
    PetscInt nStatesInput;
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

  this->data_.restoreParameterValues();

  VLOG(1) << "at end of cellml_adapter, intermediates: " << this->data_.intermediates() << " " << *this->data_.intermediates();
  this->internalTimeStepNo_++;
}


template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
void CellmlAdapter<nStates_,nIntermediates_,FunctionSpaceType>::
prepareForGetOutputConnectorData()
{
  // make representation of intermediates global, such that field variables in outputConnectorData that share the Petsc Vec's with
  // intermediates have the correct data assigned
  LOG(DEBUG) << "Transform intermediates field variable " << this->data_.intermediates() << " to global representation in order to transfer them to other solver.";
  VLOG(1) << *this->data_.intermediates();
  this->data_.intermediates()->setRepresentationGlobal();
  VLOG(1) << *this->data_.intermediates();
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

template<int nStates_, int nIntermediates_, typename FunctionSpaceType>
void CellmlAdapter<nStates_,nIntermediates_,FunctionSpaceType>::
initializeToEquilibriumValues(std::array<double,nStates_> &statesInitialValues)
{
  if (!this->rhsRoutineSingleInstance_)
  {
    LOG(ERROR) << "rhsRoutineSingleInstance is not compiled, not initializing equilibrium values for states.";
  }

  LOG(DEBUG) << "initializeToEquilibriumValues";
  LOG(INFO) << "Computing equilibrium values of states for model \"" << this->cellmlSourceCodeGenerator_.sourceFilename() << "\"...";

  double currentTime = 0.0;
  double maximumIncrement = 0;
  double dt = this->initializeStatesToEquilibriumTimestepWidth_;
  std::array<double,nStates_> previousU = statesInitialValues;
  std::array<double,nStates_> &u = statesInitialValues;
  std::array<double,nIntermediates_> intermediates;

  std::array<double,nStates_> k1;
  std::array<double,nStates_> u2, k2;
  std::array<double,nStates_> u3, k3;
  std::array<double,nStates_> u4, k4;

  int maxRateNo = 0;
  const int nInterations = 1e7;
  for (int iterationNo = 0; iterationNo < nInterations; iterationNo++)
  {
    // compute k1 = f(t, u)
    this->rhsRoutineSingleInstance_((void *)this, currentTime, u.data(), k1.data(), intermediates.data(), this->data_.parameterValues());

    // compute k2 = f(t+dt/2, u+dt/2.*k1)
    for (int stateNo = 0; stateNo < nStates_; stateNo++)
      u2[stateNo] = u[stateNo] + dt/2. * k1[stateNo];
    this->rhsRoutineSingleInstance_((void *)this, currentTime+dt/2., u2.data(), k2.data(), intermediates.data(), this->data_.parameterValues());

    // compute k3 = f(t+dt/2, u+dt/2.*k2)
    for (int stateNo = 0; stateNo < nStates_; stateNo++)
      u3[stateNo] = u[stateNo] + dt/2. * k2[stateNo];
    this->rhsRoutineSingleInstance_((void *)this, currentTime+dt/2., u3.data(), k3.data(), intermediates.data(), this->data_.parameterValues());

    // compute k4 = f(t+dt, u+dt*k3)
    for (int stateNo = 0; stateNo < nStates_; stateNo++)
      u4[stateNo] = u[stateNo] + dt * k3[stateNo];
    this->rhsRoutineSingleInstance_((void *)this, currentTime+dt, u4.data(), k4.data(), intermediates.data(), this->data_.parameterValues());

    // compute u += dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    maximumIncrement = 0;
    for (int stateNo = 0; stateNo < nStates_; stateNo++)
    {
      double increment = (k1[stateNo] + 2*k2[stateNo] + 2*k3[stateNo] + k4[stateNo])/6.;
      if (fabs(increment) > maximumIncrement)
      {
        maximumIncrement = fabs(increment);
        maxRateNo = stateNo;
      }

      u[stateNo] += dt * increment;
    }

    // if iteration diverges, restart from initial values with half the current timestep width
    if (maximumIncrement > 1e4)
    {
      LOG(WARNING) << "Search for equilibrium of states with dt=" << dt << " diverged, restarting with dt=" << dt/2.0 << "\n"
        << "Decrease value of \"initializeStatesToEquilibriumTimestepWidth\".";
      dt /= 2.0;
      u = previousU;
    }

    if (maximumIncrement < 1e-5)
    {
      LOG(INFO) << "Determined equilibrium of states after " << iterationNo << " iterations, dt=" << dt;
      break;
    }
  }

  if (maximumIncrement >= 1e-5)
  {
    LOG(ERROR) << "Equilibrium values for states were not found within " << nInterations << " RK-4 iterations! Last increment of a state: " << maximumIncrement << ", Last timestep width: " << dt;
  }

  // write computed equilibrium values to a file
  if (this->functionSpace_->meshPartition()->rankSubset()->ownRankNo() == 0)
  {

    std::stringstream filename;
    filename << this->cellmlSourceCodeGenerator_.sourceFilename() << "_equilibrium_values.txt";
    std::ofstream file(filename.str().c_str());
    if (file.is_open())
    {
      // time stamp
      auto t = std::time(nullptr);
      auto tm = *std::localtime(&t);
      std::string timeString = StringUtility::timeToString(&tm);
      file << "// Result of computation of equilibrium values for the states by opendihu on " << timeString << "\n"
        << "// Number of iterations: " << nInterations << ", dt: " << dt << "\n"
        << "// Maximum ∂u/∂t = " << maximumIncrement << " for state " << maxRateNo << "\n"
        << "// (If this is a high value, it indicates that the equilibrium was not fully reached.)\n\n";

      for (int stateNo = 0; stateNo < nStates_; stateNo++)
      {
        double lastIncrement = (k1[stateNo] + 2*k2[stateNo] + 2*k3[stateNo] + k4[stateNo])/6.;

        std::stringstream line;
        line << "state[" << stateNo << "] = " << u[stateNo] << ";";
        file << line.str() << std::string(26-line.str().length(),' ') << "// residuum: " << lastIncrement << "\n";
      }

      file << "\n  Line to copy for settings:\n  \"statesInitialValues\": [";
      for (int stateNo = 0; stateNo < nStates_; stateNo++)
      {
        if (stateNo != 0)
          file << ", ";
        file << u[stateNo];
      }
      file << "],\n";


      file.close();
      LOG(INFO) << "Values were written to \"" << filename.str() << "\".";
    }
  }
}
