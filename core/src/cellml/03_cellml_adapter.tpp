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

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
CellmlAdapter<nStates_,nAlgebraics_,FunctionSpaceType>::
CellmlAdapter(DihuContext context) :
  CallbackHandler<nStates_,nAlgebraics_,FunctionSpaceType>(context),
  Splittable()
{
  LOG(TRACE) << "CellmlAdapter constructor";

  // initialize filename in stimulation logging class from current settings
  Control::StimulationLogging logging(this->specificSettings_);
}

//! constructor from other CellmlAdapter, preserves the outputManager and context
template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
CellmlAdapter<nStates_,nAlgebraics_,FunctionSpaceType>::
CellmlAdapter(const CellmlAdapter &rhs, std::shared_ptr<FunctionSpace> functionSpace) :
  CallbackHandler<nStates_,nAlgebraics_,FunctionSpaceType>(rhs.context_, true),
  Splittable()
{
  LOG(TRACE) << "CellmlAdapter constructor from rhs";

  this->functionSpace_ = functionSpace;
  this->outputWriterManager_ = rhs.outputWriterManager_;
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
constexpr int CellmlAdapter<nStates_,nAlgebraics_,FunctionSpaceType>::
nStates()
{
  return nStates_;
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
void CellmlAdapter<nStates_,nAlgebraics_,FunctionSpaceType>::
reset()
{
  this->internalTimeStepNo_ = 0;
}
  
template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
void CellmlAdapter<nStates_,nAlgebraics_,FunctionSpaceType>::
initialize()
{
  LOG(TRACE) << "CellmlAdapter<nStates_,nAlgebraics_,FunctionSpaceType>::initialize";

  Control::PerformanceMeasurement::start("durationInitCellml");

  CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::initialize();
  
  // initialize output writers
  this->outputWriterManager_.initialize(this->context_, this->specificSettings_, this->data().functionSpace()->meshPartition()->rankSubset());

  // load rhs routine
  this->initializeRhsRoutine();

  // parse the callback functions from the python config
  this->initializeCallbackFunctions();
  
  Control::PerformanceMeasurement::stop("durationInitCellml");

  if (this->specificSettings_.hasKey("prefactor"))
  {
    LOG(WARNING) << this->specificSettings_ << "[\"prefactor\"] is no longer an option! There is no more functionality to scale values during transfer.";
  }
  if (this->outputWriterManager_.hasOutputWriters())
  {
    LOG(INFO) << "CellML has output writers. This will be slow if it outputs lots of data and should only be used for debugging.";
  }

  this->internalTimeStepNo_ = 0;

  this->setSpecificStatesCallFrequency_ = this->specificSettings_.getOptionDouble("setSpecificStatesCallFrequency", 0.0);
  this->specificSettings_.getOptionVector("setSpecificStatesFrequencyJitter", this->setSpecificStatesFrequencyJitter_);
  this->setSpecificStatesRepeatAfterFirstCall_ = this->specificSettings_.getOptionDouble("setSpecificStatesRepeatAfterFirstCall", 0.0);
  this->setSpecificStatesCallEnableBegin_ = this->specificSettings_.getOptionDouble("setSpecificStatesCallEnableBegin", 0.0);

  // if setSpecificStatesCallFrequency_ is set to None, set to 0
  if (this->setSpecificStatesCallFrequency_ == std::numeric_limits<double>::max())
    this->setSpecificStatesCallFrequency_ = 0;

  // if setSpecificStatesRepeatAfterFirstCall_ is set to None, set to 0
  if (this->setSpecificStatesRepeatAfterFirstCall_ == std::numeric_limits<double>::max())
    this->setSpecificStatesRepeatAfterFirstCall_ = 0;

  if (this->setSpecificStatesCallFrequency_ != 0 && this->setSpecificStatesRepeatAfterFirstCall_ == 0)
  {
    LOG(FATAL) << "In " << this->specificSettings_ << ", you have set \"setSpecificStatesCallFrequency\", but not "
      << "\"setSpecificStatesRepeatAfterFirstCall\". This is not allowed.\n"
      << "Either set \"setSpecificStatesRepeatAfterFirstCall\" to a reasonable value or do not use \"setSpecificStatesCallFrequency\" at "
      << "all (set to None or 0 and instead use \"setSpecificStatesCallInterval\").";
  }

  // initialize the lastCallSpecificStatesTime_
  this->currentJitter_ = 0;
  this->jitterIndex_ = 0;
  this->lastCallSpecificStatesTime_ = this->setSpecificStatesCallEnableBegin_ - 1e-13 - 1./(this->setSpecificStatesCallFrequency_+this->currentJitter_);

  LOG(DEBUG) << "Cellml end of initialize, " << this->sourceToCompileFilename_ << ", statesForTransfer: " << this->data_.statesForTransfer();
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
void CellmlAdapter<nStates_,nAlgebraics_,FunctionSpaceType>::
initializeForImplicitTimeStepping()
{
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
void CellmlAdapter<nStates_,nAlgebraics_,FunctionSpaceType>::
setRankSubset(Partition::RankSubset rankSubset)
{
  // do nothing because we don't have stored data here (the data on which the computation is performed comes in evaluateTimesteppingRightHandSide from parameters) 
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
void CellmlAdapter<nStates_,nAlgebraics_,FunctionSpaceType>::
evaluateTimesteppingRightHandSideExplicit(Vec& input, Vec& output, int timeStepNo, double currentTime)
{
  // prepare variable to work on
  // get raw pointers from Petsc data structures
  double *statesLocal;
  double *ratesLocal;
  double *algebraicsLocal;
  PetscErrorCode ierr;
  ierr = VecGetArray(input, &statesLocal); CHKERRV(ierr);   // get r/w pointer to contiguous array of the data, VecRestoreArray() needs to be called afterwards
  ierr = VecGetArray(output, &ratesLocal); CHKERRV(ierr);
  ierr = VecGetArray(this->data_.algebraics()->getValuesContiguous(), &algebraicsLocal); CHKERRV(ierr);

  // get sizes of input and output Vecs
  PetscInt nStatesInput, nRates, nAlgebraics = 101;
  ierr = VecGetSize(input, &nStatesInput); CHKERRV(ierr);
  ierr = VecGetSize(output, &nRates); CHKERRV(ierr);
  ierr = VecGetLocalSize(this->data_.algebraics()->getValuesContiguous(), &nAlgebraics); CHKERRV(ierr);

  LOG(DEBUG) << "Cellml evaluateTimesteppingRightHandSideExplicit, " << this->sourceToCompileFilename_ << ", statesForTransfer: " << this->data_.statesForTransfer();
  nAlgebraics = nAlgebraics/this->nInstances_;
  VLOG(1) << "algebraics array has " << nAlgebraics << " entries";

  // check validity of sizes
  VLOG(1) << "evaluateTimesteppingRightHandSideExplicit, input nStates_: " << nStatesInput << ", output nRates: " << nRates;
  VLOG(1) << "timeStepNo: " << timeStepNo << ", currentTime: " << currentTime << ", internalTimeStepNo: " << this->internalTimeStepNo_;

  if (nStatesInput != nStates_*this->nInstances_)
  {
    LOG(ERROR) << "nStatesInput does not match nStates and nInstances! nStatesInput=" << nStatesInput << ", nStates_=" << nStates_ << ", nInstances=" << this->nInstances_;
  }
  assert (nStatesInput == nStates_*this->nInstances_);
  assert (nRates == nStates_*this->nInstances_);

  if (nAlgebraics != nAlgebraics_)
  {
    LOG(FATAL) << "nInstances: " << this->nInstances_ << ", nAlgebraics (size of vector / nInstances): " << nAlgebraics << ", nAlgebraics_: " << nAlgebraics_;
  }
  assert (nAlgebraics == nAlgebraics_);

  // make the parameterValues_ vector available
  this->data_.prepareParameterValues();

  //LOG(DEBUG) << " evaluateTimesteppingRightHandSide: nInstances=" << this->nInstances_ << ", nStates_=" << nStates_;

  // handle callback functions "setSpecificParameters" and "setSpecificStates"
  checkCallbackParameters(currentTime);
  checkCallbackStates(currentTime, statesLocal);

#if 0
  if (this->sourceToCompileFilename_.substr(0, 11) == "src/spindle")
  {
    CLOG(INFO, "special") << "CellML \"" << this->sourceToCompileFilename_ << "\" parameters:";
    for (int instanceNo = 0; instanceNo < this->nInstances_; instanceNo++)
    {
      for (int parameterNo = 0; parameterNo < this->cellmlSourceCodeGenerator_.nParameters(); parameterNo++)
      {
        CLOG(INFO, "special") << "  instance " << instanceNo << " parameter " << parameterNo << ": " << this->data_.parameterValues()[parameterNo*this->nInstances_ + instanceNo];
      }
    }
  }
#endif

  // call actual rhs method
  if (this->rhsRoutine_)
  {
    //Control::PerformanceMeasurement::start("rhsEvaluationTime");  // commented out because it takes too long in this very inner loop

    // call actual rhs routine from cellml code
    this->rhsRoutine_((void *)this, currentTime, statesLocal, ratesLocal, algebraicsLocal, this->data_.parameterValues());

    //Control::PerformanceMeasurement::stop("rhsEvaluationTime");
  }

  // handle callback function "handleResult"
  checkCallbackAlgebraics(currentTime, statesLocal, algebraicsLocal);

  // give control of data back to Petsc
  ierr = VecRestoreArray(input, &statesLocal); CHKERRV(ierr);
  ierr = VecRestoreArray(output, &ratesLocal); CHKERRV(ierr);
  ierr = VecRestoreArray(this->data_.algebraics()->getValuesContiguous(), &algebraicsLocal); CHKERRV(ierr);

  this->data_.restoreParameterValues();

  LOG(DEBUG) << "Cellml end of evaluateTimesteppingRightHandSideExplicit, " << this->sourceToCompileFilename_ << ", statesForTransfer: " << this->data_.statesForTransfer();
  // call output writer to write output files
  this->outputWriterManager_.writeOutput(this->data_, this->internalTimeStepNo_, currentTime);

  VLOG(1) << "at end of cellml_adapter, algebraics: " << this->data_.algebraics() << " " << *this->data_.algebraics();
  this->internalTimeStepNo_++;
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
void CellmlAdapter<nStates_,nAlgebraics_,FunctionSpaceType>::
checkCallbackParameters(double currentTime)
{
  // get new values for parameters, call callback function of python config
  if (this->pythonSetSpecificParametersFunction_ && this->internalTimeStepNo_ % this->setSpecificParametersCallInterval_ == 0)
  {
    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;

    VLOG(1) << "call setSpecificParameters";
    this->callPythonSetSpecificParametersFunction(this->nInstances_, this->internalTimeStepNo_, currentTime, this->data_.parameterValues(), this->cellmlSourceCodeGenerator_.nParameters());
  }
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
void CellmlAdapter<nStates_,nAlgebraics_,FunctionSpaceType>::
checkCallbackStates(double currentTime, double *statesLocal)
{
  // there is a parallel piece of code to this one in FastMonodomainSolverBase<>::isCurrentPointStimulated(), specialized_solver/fast_monodomain_solver/fast_monodomain_solver_compute.tpp

  VLOG(1) << "currentTime: " << currentTime << ", lastCallSpecificStatesTime_: " << this->lastCallSpecificStatesTime_
    << ", setSpecificStatesCallFrequency_: " << this->setSpecificStatesCallFrequency_ << ", "
    << this->lastCallSpecificStatesTime_ + 1./this->setSpecificStatesCallFrequency_;
  VLOG(1) << "this->pythonSetSpecificStatesFunction_? " << (this->pythonSetSpecificStatesFunction_? "true" : "false")
    << ", this->setSpecificStatesCallInterval_: " << this->setSpecificStatesCallInterval_ << " != 0? " << (this->setSpecificStatesCallInterval_ != 0? "true" : "false")
    << ", (this->internalTimeStepNo_ % this->setSpecificStatesCallInterval_) = " << this->internalTimeStepNo_  << " % " << this->setSpecificStatesCallInterval_
    << ", this->setSpecificStatesCallFrequency_= " << this->setSpecificStatesCallFrequency_ << " != 0.0? " << (this->setSpecificStatesCallFrequency_ != 0.0? "true" : "false")
    << ", currentTime=" << currentTime << " >= " << this->lastCallSpecificStatesTime_ << " + " << 1./this->setSpecificStatesCallFrequency_ << " = " << this->lastCallSpecificStatesTime_ + 1./this->setSpecificStatesCallFrequency_ << "? "
    << (currentTime >= this->lastCallSpecificStatesTime_ + 1./this->setSpecificStatesCallFrequency_? "true" : "false");

  bool stimulate = false;

  // get new values for states, call callback function of python config
  if (this->pythonSetSpecificStatesFunction_
      && (
          (this->setSpecificStatesCallInterval_ != 0 && this->internalTimeStepNo_ % this->setSpecificStatesCallInterval_ == 0)
          || (this->setSpecificStatesCallFrequency_ != 0.0 && currentTime >= this->lastCallSpecificStatesTime_ + 1./(this->setSpecificStatesCallFrequency_+this->currentJitter_)
              && currentTime >= this->setSpecificStatesCallEnableBegin_-1e-13)
         )
     )
  {
    stimulate = true;

    // if current stimulation is over
    if (this->setSpecificStatesRepeatAfterFirstCall_ != 0
        && currentTime - (this->lastCallSpecificStatesTime_ + 1./(this->setSpecificStatesCallFrequency_+this->currentJitter_)) > this->setSpecificStatesRepeatAfterFirstCall_)
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

  VLOG(1) << "stimulate = " << stimulate;

  static bool currentlyStimulating = false;
  if (stimulate)
  {
    VLOG(1) << "currentlyStimulating: " << currentlyStimulating;

    // if this is the first point in time of the current stimulation, log stimulation time
    if (!currentlyStimulating)
    {
      currentlyStimulating = true;
      Control::StimulationLogging::logStimulationBegin(currentTime, -1, this->fiberNoGlobal_);
    }

    VLOG(1) << "call setSpecificStates, this->internalTimeStepNo_ = " << this->internalTimeStepNo_ << ", this->setSpecificStatesCallInterval_: " << this->setSpecificStatesCallInterval_;
    VLOG(1) << "currentTime: " << currentTime << ", call setSpecificStates, this->internalTimeStepNo_ = " << this->internalTimeStepNo_ << ", this->setSpecificStatesCallInterval_: " << this->setSpecificStatesCallInterval_;
    this->callPythonSetSpecificStatesFunction(this->nInstances_, this->internalTimeStepNo_, currentTime, statesLocal);
  }
  else
  {
    currentlyStimulating = false;
  }
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
void CellmlAdapter<nStates_,nAlgebraics_,FunctionSpaceType>::
checkCallbackAlgebraics(double currentTime, double *statesLocal, double *algebraicsLocal)
{
  // handle resulting algebraics, call callback function of python config
  if (this->pythonHandleResultFunction_ && this->internalTimeStepNo_ % this->handleResultCallInterval_ == 0)
  {
    // start critical section for python API calls
    // PythonUtility::GlobalInterpreterLock lock;

    this->callPythonHandleResultFunction(this->nInstances_, this->internalTimeStepNo_, currentTime, statesLocal, algebraicsLocal);
  }
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
void CellmlAdapter<nStates_,nAlgebraics_,FunctionSpaceType>::
prepareForGetOutputConnectorData()
{
  // This method is called before getOutputConnectorData() of the timestepping scheme.

  // make representation of algebraics global, such that field variables in outputConnectorData that share the Petsc Vec's with
  // algebraics have the correct data assigned
  LOG(DEBUG) << "Transform algebraics and parameters field variables to global representation in order to transfer them to other solver, such that extracted component-field variables in timestepping scheme have the correct values.";

  VLOG(1) << *this->data_.algebraics();
  this->data_.algebraics()->setRepresentationGlobal();
  VLOG(1) << *this->data_.algebraics();

  VLOG(1) << *this->data_.parameters();
  this->data_.parameters()->setRepresentationGlobal();
  VLOG(1) << *this->data_.parameters();
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
void CellmlAdapter<nStates_,nAlgebraics_,FunctionSpaceType>::
getComponentNames(std::vector<std::string> &stateNames)
{
  this->getStateNames(stateNames);
}

//! return the mesh
template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
std::shared_ptr<FunctionSpaceType> CellmlAdapter<nStates_,nAlgebraics_,FunctionSpaceType>::
functionSpace()
{
  return CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::functionSpace();
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
template<typename FunctionSpaceType2>
bool CellmlAdapter<nStates_,nAlgebraics_,FunctionSpaceType>::
setInitialValues(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nStates_>> initialValues)
{
  return CellmlAdapterBase<nStates_,nAlgebraics_,FunctionSpaceType>::template setInitialValues<FunctionSpaceType2>(initialValues);
}

template<int nStates_, int nAlgebraics_, typename FunctionSpaceType>
void CellmlAdapter<nStates_,nAlgebraics_,FunctionSpaceType>::
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
  std::array<double,nAlgebraics_> algebraics;

  std::array<double,nStates_> k1;
  std::array<double,nStates_> u2, k2;
  std::array<double,nStates_> u3, k3;
  std::array<double,nStates_> u4, k4;

  int maxRateNo = 0;
  const int nInterations = 1e7;
  for (int iterationNo = 0; iterationNo < nInterations; iterationNo++)
  {
    // compute k1 = f(t, u)
    this->rhsRoutineSingleInstance_((void *)this, currentTime, u.data(), k1.data(), algebraics.data(), this->data_.parameterValues());

    // compute k2 = f(t+dt/2, u+dt/2.*k1)
    for (int stateNo = 0; stateNo < nStates_; stateNo++)
      u2[stateNo] = u[stateNo] + dt/2. * k1[stateNo];
    this->rhsRoutineSingleInstance_((void *)this, currentTime+dt/2., u2.data(), k2.data(), algebraics.data(), this->data_.parameterValues());

    // compute k3 = f(t+dt/2, u+dt/2.*k2)
    for (int stateNo = 0; stateNo < nStates_; stateNo++)
      u3[stateNo] = u[stateNo] + dt/2. * k2[stateNo];
    this->rhsRoutineSingleInstance_((void *)this, currentTime+dt/2., u3.data(), k3.data(), algebraics.data(), this->data_.parameterValues());

    // compute k4 = f(t+dt, u+dt*k3)
    for (int stateNo = 0; stateNo < nStates_; stateNo++)
      u4[stateNo] = u[stateNo] + dt * k3[stateNo];
    this->rhsRoutineSingleInstance_((void *)this, currentTime+dt, u4.data(), k4.data(), algebraics.data(), this->data_.parameterValues());

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
