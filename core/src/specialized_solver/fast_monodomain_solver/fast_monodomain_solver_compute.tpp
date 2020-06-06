#include "specialized_solver/fast_monodomain_solver/fast_monodomain_solver_base.h"

#include "partition/rank_subset.h"
#include "control/diagnostic_tool/stimulation_logging.h"

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
run()
{
  initialize();
  advanceTimeSpan();
}

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
advanceTimeSpan()
{
  LOG(TRACE) << "FastMonodomainSolver::advanceTimeSpan";

  // loop over fibers and communicate element lengths and initial values to the ranks that participate in computing
  fetchFiberData();

  //Control::PerformanceMeasurement::startFlops();

  // do computation of own fibers, HH is hardcoded, stimulation from parsed MU and firing_times files
  computeMonodomain();

  //Control::PerformanceMeasurement::endFlops();

  // loop over fibers and communicate resulting values back
  updateFiberData();

  //nestedSolvers_.advanceTimeSpan();

  // call output writer of diffusion
  std::vector<typename NestedSolversType::TimeSteppingSchemeType> &instances = nestedSolvers_.instancesLocal();

  for (int i = 0; i < instances.size(); i++)
  {
    // call write output of MultipleInstances, callCountIncrement is the number of times the output writer would have been called without FastMonodomainSolver
    instances[i].timeStepping2().writeOutput(0, currentTime_, nTimeStepsSplitting_);
  }

  // call own output writers, write current 0D output values, not yet implemented
  //this->outputWriterManager_.writeOutput(*this->data_, 0, currentTime_);
}

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
computeMonodomain()
{
  LOG(TRACE) << "computeMonodomain";

  // initialize data vector
  // array of vectorized struct

  // fetch timestep widths and total time span
  std::vector<typename NestedSolversType::TimeSteppingSchemeType> &instances = nestedSolvers_.instancesLocal();

  TimeSteppingScheme::Heun<CellmlAdapterType> &heun = instances[0].timeStepping1().instancesLocal()[0];
  durationLogKey0D_ = heun.durationLogKey();

  DiffusionTimeSteppingScheme &implicitEuler = instances[0].timeStepping2().instancesLocal()[0];
  durationLogKey1D_ = implicitEuler.durationLogKey();
  double prefactor = implicitEuler.discretizableInTime().data().context().getPythonConfig().getOptionDouble("prefactor", 1.0);

  LOG(DEBUG) << "durationLogKeys: " << durationLogKey0D_ << "," << durationLogKey1D_;

  double startTime = instances[0].startTime();
  double timeStepWidthSplitting = instances[0].timeStepWidth();
  nTimeStepsSplitting_ = instances[0].numberTimeSteps();

  heun.setTimeSpan(startTime, startTime + 0.5 * timeStepWidthSplitting);
  double dt0D = heun.timeStepWidth();
  int nTimeSteps0D = heun.numberTimeSteps();

  implicitEuler.setTimeSpan(startTime, startTime + timeStepWidthSplitting);
  double dt1D = implicitEuler.timeStepWidth();
  int nTimeSteps1D = implicitEuler.numberTimeSteps();

  LOG(DEBUG) << "prefactor: " << prefactor << ", dtSplitting: " << timeStepWidthSplitting << ", n steps: " << nTimeStepsSplitting_;
  LOG(DEBUG) << "dt0D: " << dt0D << ", n steps: " << nTimeSteps0D << ", dt1D: " << dt1D << ", n steps: " << nTimeSteps1D;

  // picture for strang splitting:
  // ===== t ==>
  // -1->
  //  /
  // ---2--->
  //      /
  //     -1->
  //        |
  //        2

  // loop over splitting time steps
  for (int timeStepNo = 0; timeStepNo < nTimeStepsSplitting_; timeStepNo++)
  {
    // perform Strang splitting
    double currentTime = startTime + timeStepNo * timeStepWidthSplitting;

    LOG(DEBUG) << "splitting " << timeStepNo << "/" << nTimeStepsSplitting_ << ", t: " << currentTime;

    // compute midTime once per step to reuse it. [currentTime, midTime=currentTime+0.5*timeStepWidth, currentTime+timeStepWidth]
    double midTime = currentTime + 0.5 * timeStepWidthSplitting;
    bool storeAlgebraicsForTransfer = timeStepNo == nTimeStepsSplitting_-1;   // after the last timestep, store the algebraics for transfer

    // perform splitting
    compute0D(currentTime, dt0D, nTimeSteps0D, false);
    compute1D(currentTime, dt1D, nTimeSteps1D, prefactor);
    compute0D(midTime,     dt0D, nTimeSteps0D, storeAlgebraicsForTransfer);
  }

  currentTime_ = instances[0].endTime();
}

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
compute0D(double startTime, double timeStepWidth, int nTimeSteps, bool storeAlgebraicsForTransfer)
{
  Control::PerformanceMeasurement::start(durationLogKey0D_);
  LOG(DEBUG) << "compute0D(" << startTime << "), " << nTimeSteps << " time step" << (nTimeSteps == 1? "" : "s");

  using Vc::double_v;

  // Heun scheme:
  // y* = y_n + dt*rhs(y_n)
  // y_n+1 = y_n + 0.5*[rhs(y_n) + rhs(y*)]

  // loop over point buffers, i.e., sets of 4 neighouring points of the fiber
  int nPointBuffers = fiberPointBuffers_.size();
  for (global_no_t pointBuffersNo = 0; pointBuffersNo < nPointBuffers; pointBuffersNo++)
  {
    int fiberDataNo = pointBuffersNo * Vc::double_v::Size / fiberData_[0].valuesLength;
    int indexInFiber = pointBuffersNo * Vc::double_v::Size - fiberData_[fiberDataNo].valuesOffset;

    // determine if current point is at center of fiber
    int fiberCenterIndex = fiberData_[fiberDataNo].valuesLength / 2;
    bool currentPointIsInCenter = ((fiberCenterIndex - indexInFiber) < Vc::double_v::Size);
    VLOG(3) << "currentPointIsInCenter: " << currentPointIsInCenter << ", pointBuffersNo: " << pointBuffersNo << ", fiberDataNo: " << fiberDataNo << ", indexInFiber:" << indexInFiber << ", fiberCenterIndex: " << fiberCenterIndex << ", " << (indexInFiber - fiberCenterIndex) << " < " << Vc::double_v::Size;

    VLOG(3) << "pointBuffersNo: " << pointBuffersNo << ", fiberDataNo: " << fiberDataNo << ", indexInFiber: " << indexInFiber;

    // save previous state values for equilibrium acceleration
    Vc::double_v statesPreviousValues[nStates];

    if (disableComputationWhenStatesAreCloseToEquilibrium_)
    {
      for (int stateNo = 0; stateNo < nStates; stateNo++)
      {
        statesPreviousValues[stateNo] = fiberPointBuffers_[pointBuffersNo].states[stateNo];
      }
    }

    // loop over timesteps
    for (int timeStepNo = 0; timeStepNo < nTimeSteps; timeStepNo++)
    {
      // determine if fiber gets stimulated
      double currentTime = startTime + timeStepNo * timeStepWidth;

      // check if current point will be stimulated
      const bool stimulateCurrentPoint = isCurrentPointStimulated(fiberDataNo, currentTime, currentPointIsInCenter);
      const bool argumentStoreAlgebraics = storeAlgebraicsForTransfer && timeStepNo == nTimeSteps-1;

      if (stimulateCurrentPoint)
      {
        LOG(INFO) << "t: " << currentTime << ", stimulate fiber " << fiberData_[fiberDataNo].fiberNoGlobal
          << ", MU " << fiberData_[fiberDataNo].motorUnitNo;
      }

      // if the current point does not need to get computed because the value won't change
      if (isEquilibriumAccelerationCurrentPointDisabled(stimulateCurrentPoint, pointBuffersNo))
      {
        continue;
      }

      // do not compute fiber if respective option is set and the fiber has not yet been stimulated
      if (onlyComputeIfHasBeenStimulated_ && !fiberHasBeenStimulated_[fiberDataNo])
      {
        continue;
      }

      // call method to compute 0D problem
      assert (compute0DInstance_ != nullptr);
      compute0DInstance_(fiberPointBuffers_[pointBuffersNo].states, fiberPointBuffersParameters_[pointBuffersNo],
                         currentTime, timeStepWidth, stimulateCurrentPoint,
                         argumentStoreAlgebraics, fiberPointBuffersAlgebraicsForTransfer_[pointBuffersNo],
                         algebraicsForTransfer_, valueForStimulatedPoint_);
    }  // loop over timesteps

    equilibriumAccelerationUpdate(statesPreviousValues, pointBuffersNo);

    //VLOG(3) << "-> index " << pointBuffersNo << ", states [" << state0 << "," << state1 << "," << state2 << "," << state3 << "]";
  }

  // visualize equilibrium states for debugging
#if 0
  if (storeAlgebraicsForTransfer)
  {
    std::stringstream s;
    int fiberDataNoPrevious = -1;
    for (global_no_t pointBuffersNo = 0; pointBuffersNo < nPointBuffers; pointBuffersNo++)
    {
      int fiberDataNo = pointBuffersNo * Vc::double_v::Size / fiberData_[0].valuesLength;
      if (fiberDataNoPrevious != fiberDataNo)
      {
        s << std::endl << fiberDataNo;
        fiberDataNoPrevious = fiberDataNo;
      }

      switch(fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo])
      {
      case constant:
        s << ".";
        break;
      case neighbour_not_constant:
        s << "n";
        break;
      case not_constant:
        s << "N";
        break;
      }
    }
    LOG(INFO) << s.str();
  }
#endif

  VLOG(1) << "nFiberPointBuffers: " << nPointBuffers;
  Control::PerformanceMeasurement::stop(durationLogKey0D_);
}

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
compute1D(double startTime, double timeStepWidth, int nTimeSteps, double prefactor)
{
  // if all entries are at equilibrium, nothing will be computed, skip also 1D computation
  if (nFiberPointBufferStatesCloseToEquilibrium_ == fiberPointBuffersStatesAreCloseToEquilibrium_.size())
  {
    LOG(DEBUG) << "skip 1D computation";
    return;
  }

  Control::PerformanceMeasurement::start(durationLogKey1D_);

  LOG(DEBUG) << "compute1D(" << startTime << ")";

  // depending on DiffusionTimeSteppingScheme either do Implicit Euler or Crank-Nicolson
  // Implicit Euler step:
  // (K - 1/dt*M) u^{n+1} = -1/dt*M u^{n})
  // Crank-Nicolson step:
  // (1/2*K - 1/dt*M) u^{n+1} = (-1/2*K -1/dt*M) u^{n})

  // stencil K: 1/h*[_-1_  1  ]*prefactor
  // stencil M:   h*[_1/3_ 1/6]

  bool useImplicitEuler = std::is_same<DiffusionTimeSteppingScheme,
                            TimeSteppingScheme::ImplicitEuler<typename DiffusionTimeSteppingScheme::DiscretizableInTime>
                          >::value;

  const double dt = timeStepWidth;

  // loop over fibers
  for (int fiberDataNo = 0; fiberDataNo < fiberData_.size(); fiberDataNo++)
  {
    int nValues = fiberData_[fiberDataNo].vmValues.size();

#ifndef NDEBUG
    VLOG(1) << "fiber " << fiberDataNo << "/" << fiberData_.size() << ", valuesOffset: " << fiberData_[fiberDataNo].valuesOffset
      << ", has " << nValues << " values: ";

    std::stringstream s, s2;
    for (int valueNo = 0; valueNo < nValues; valueNo++)
    {
      global_no_t valuesIndexAllFibers = fiberData_[fiberDataNo].valuesOffset + valueNo;
      global_no_t pointBuffersNo = valuesIndexAllFibers / Vc::double_v::Size;
      int entryNo = valuesIndexAllFibers % Vc::double_v::Size;
      double u = fiberPointBuffers_[pointBuffersNo].states[0][entryNo];
      if (valueNo != 0)
      {
        s << ", ";
        s2 << ",";
      }
      s << u;
      s2 << "[" << pointBuffersNo << "][" << entryNo << "](" << valuesIndexAllFibers << ")  ";
    }
    VLOG(1) << s.str();
    //LOG(DEBUG) << "indices: " << s2.str();
#endif
    // [ b c     ] [x]   [d]
    // [ a b c   ] [x] = [d]
    // [   a b c ] [x]   [d]
    // [     a b ] [x]   [d]

    // Thomas algorithm
    // forward substitution
    // c'_0 = c_0 / b_0
    // c'_i = c_i / (b_i - c'_{i-1}*a_i)

    // d'_0 = d_0 / b_0
    // d'_i = (d_i - d'_{i-1}*a_i) / (b_i - c'_{i-1}*a_i)

    // backward substitution
    // x_n = d'_n
    // x_i = d'_i - c'_i * x_{i+1}

    // helper buffers c', d'
    static std::vector<double> cAlgebraic(nValues-1);
    static std::vector<double> dAlgebraic(nValues);

    // perform forward substitution
    // loop over entries / rows of matrices
    for (int valueNo = 0; valueNo < nValues; valueNo++)
    {
     // new with CN
      double a = 0;
      double b = 0;
      double c = 0;
      double d = 0;

      double u_previous = 0;
      double u_center = 0;
      double u_next = 0;

      global_no_t valuesIndexAllFibers = fiberData_[fiberDataNo].valuesOffset + valueNo;
      global_no_t pointBuffersNo = valuesIndexAllFibers / Vc::double_v::Size;
      int entryNo = valuesIndexAllFibers % Vc::double_v::Size;
      u_center = fiberPointBuffers_[pointBuffersNo].states[0][entryNo];

      // contribution from left element
      if (valueNo > 0)
      {
        global_no_t valuesIndexAllFibers = fiberData_[fiberDataNo].valuesOffset + valueNo - 1;
        global_no_t pointBuffersNo = valuesIndexAllFibers / Vc::double_v::Size;
        int entryNo = valuesIndexAllFibers % Vc::double_v::Size;
        u_previous = fiberPointBuffers_[pointBuffersNo].states[0][entryNo];

        // stencil K: 1/h*[1   _-1_ ]*prefactor
        // stencil M:   h*[1/6 _1/3_]

        double h_left = fiberData_[fiberDataNo].elementLengths[valueNo-1];
        double k_left = 1./h_left*(1) * prefactor;
        double m_left = h_left*1./6;

        if (useImplicitEuler)
        {
          a = (k_left - 1/dt*m_left);
        }
        else  // Crank-Nicolson
        {
          a = (k_left/2. - 1/dt*m_left);
        }

        double k_right = 1./h_left*(-1) * prefactor;
        double m_right = h_left*1./3;

        if (useImplicitEuler)
        {
          b += (k_right - 1/dt*m_right);
          d += (-1/dt*m_left) * u_previous + (-1/dt*m_right) * u_center;
        }
        else  // Crank-Nicolson
        {
          b += (k_right/2. - 1/dt*m_right);
          d += (-k_left/2. - 1/dt*m_left) * u_previous + (-k_right/2. - 1/dt*m_right) * u_center;
        }
      }

      // contribution from right element
      if (valueNo < nValues-1)
      {
        global_no_t valuesIndexAllFibers = fiberData_[fiberDataNo].valuesOffset + valueNo + 1;
        global_no_t pointBuffersNo = valuesIndexAllFibers / Vc::double_v::Size;
        int entryNo = valuesIndexAllFibers % Vc::double_v::Size;
        u_next = fiberPointBuffers_[pointBuffersNo].states[0][entryNo];

        // stencil K: 1/h*[_-1_  1  ]*prefactor
        // stencil M:   h*[_1/3_ 1/6]

        double h_right = fiberData_[fiberDataNo].elementLengths[valueNo];
        double k_right = 1./h_right*(1) * prefactor;
        double m_right = h_right*1./6;

        if (useImplicitEuler)
        {
          c = (k_right - 1/dt*m_right);
        }
        else  // Crank-Nicolson
        {
          c = (k_right/2. - 1/dt*m_right);
        }

        double k_left = 1./h_right*(-1) * prefactor;
        double m_left = h_right*1./3;

        if (useImplicitEuler)
        {
          b += (k_left - 1/dt*m_left);
          d += (-1/dt*m_left) * u_center + (-1/dt*m_right) * u_next;
        }
        else  // Crank-Nicolson
        {
          b += (k_left/2. - 1/dt*m_left);
          d += (-k_left/2. - 1/dt*m_left) * u_center + (-k_right/2. - 1/dt*m_right) * u_next;
        }
      }

      if (valueNo == 0)
      {
        // c'_0 = c_0 / b_0
        cAlgebraic[valueNo] = c / b;

        // d'_0 = d_0 / b_0
        dAlgebraic[valueNo] = d / b;
      }
      else
      {
        if (valueNo != nValues-1)
        {
          // c'_i = c_i / (b_i - c'_{i-1}*a_i)
          cAlgebraic[valueNo] = c / (b - cAlgebraic[valueNo-1]*a);
        }

        // d'_i = (d_i - d'_{i-1}*a_i) / (b_i - c'_{i-1}*a_i)
        dAlgebraic[valueNo] = (d - dAlgebraic[valueNo-1]*a) / (b - cAlgebraic[valueNo-1]*a);
      }

#ifndef NDEBUG
      if (valueNo < 5 || valueNo >= nValues-6)
      {
        VLOG(2) << "valueNo: " << valueNo << ", a: " << a << ", b: " << b << ", c: " << c << ", d: " << d
          << ", u: " << u_center << ", c': " << (valueNo < nValues-1? cAlgebraic[valueNo] : 0) << ", d': " << dAlgebraic[valueNo]
          << " = (" << d << "-" << dAlgebraic[valueNo-1]*a << ")/(" << b << "-" << cAlgebraic[valueNo-1]*a << ") = " << (d - dAlgebraic[valueNo-1]*a) << "/" << (b - cAlgebraic[valueNo-1]*a);
      }
#endif
    }

    //LOG(DEBUG) << "cAlgebraic: " << cAlgebraic;
    //LOG(DEBUG) << "dAlgebraic: " << dAlgebraic;

    // perform backward substitution
    // x_n = d'_n
    global_no_t valuesIndexAllFibers = fiberData_[fiberDataNo].valuesOffset + nValues-1;
    global_no_t pointBuffersNo = valuesIndexAllFibers / Vc::double_v::Size;
    int entryNo = valuesIndexAllFibers % Vc::double_v::Size;
    fiberPointBuffers_[pointBuffersNo].states[0][entryNo] = dAlgebraic[nValues-1];
    //LOG(DEBUG) << "set entry (" << pointBuffersNo << "," << entryNo << ") to " << dAlgebraic[nValues-1];

    double previousValue = dAlgebraic[nValues-1];

    // loop over entries / rows of matrices
    for (int valueNo = nValues-2; valueNo >= 0; valueNo--)
    {
      // x_i = d'_i - c'_i * x_{i+1}
      global_no_t valuesIndexAllFibers = fiberData_[fiberDataNo].valuesOffset + valueNo;
      global_no_t pointBuffersNo = valuesIndexAllFibers / Vc::double_v::Size;
      int entryNo = valuesIndexAllFibers % Vc::double_v::Size;

      double resultValue = dAlgebraic[valueNo] - cAlgebraic[valueNo] * previousValue;
      fiberPointBuffers_[pointBuffersNo].states[0][entryNo] = resultValue;

      previousValue = resultValue;
    }

    int nPointBuffers = fiberPointBuffers_.size();
    pointBuffersNo = int(nPointBuffers/2);

#ifndef NDEBUG
    s.str("");
    for (int valueNo = 0; valueNo < nValues; valueNo++)
    {
      global_no_t valuesIndexAllFibers = fiberData_[fiberDataNo].valuesOffset + valueNo;
      global_no_t pointBuffersNo = valuesIndexAllFibers / Vc::double_v::Size;
      int entryNo = valuesIndexAllFibers % Vc::double_v::Size;
      double u = fiberPointBuffers_[pointBuffersNo].states[0][entryNo];
      if (valueNo != 0)
        s << ", ";
      s << u;
    }
    VLOG(1) << " -> " << s.str();
#endif
  }
  Control::PerformanceMeasurement::stop(durationLogKey1D_);
}

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
bool FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
isCurrentPointStimulated(int fiberDataNo, double currentTime, bool currentPointIsInCenter)
{
  FiberData &fiberDataCurrentPoint = fiberData_[fiberDataNo];

  // there is a parallel piece of code to this one in CellmlAdapter<>::checkCallbackStates, cellml/03_cellml_adapter.tpp

  // get time from with testing for stimulation is enabled
  const double lastStimulationCheckTime                 = fiberDataCurrentPoint.lastStimulationCheckTime;

  const double setSpecificStatesCallFrequency           = fiberDataCurrentPoint.setSpecificStatesCallFrequency;
  const double setSpecificStatesRepeatAfterFirstCall    = fiberDataCurrentPoint.setSpecificStatesRepeatAfterFirstCall;
  const double setSpecificStatesCallEnableBegin         = fiberDataCurrentPoint.setSpecificStatesCallEnableBegin;

  const int &motorUnitNo                                = fiberDataCurrentPoint.motorUnitNo;
  std::vector<double> &setSpecificStatesFrequencyJitter = fiberDataCurrentPoint.setSpecificStatesFrequencyJitter;
  int &jitterIndex                                      = fiberDataCurrentPoint.jitterIndex;
  double &currentJitter                                 = fiberDataCurrentPoint.currentJitter;

  // check if time has come to call setSpecificStates
  bool checkStimulation = false;

  VLOG(1) << "currentTime: " << currentTime << ", lastStimulationCheckTime: " << lastStimulationCheckTime << ", next time point: " << lastStimulationCheckTime + 1./(setSpecificStatesCallFrequency+currentJitter);
  VLOG(1) << "setSpecificStatesCallFrequency: " << setSpecificStatesCallFrequency << ", currentJitter: " << currentJitter << ", setSpecificStatesCallEnableBegin: " << setSpecificStatesCallEnableBegin;

  if (currentTime >= lastStimulationCheckTime + 1./(setSpecificStatesCallFrequency+currentJitter)
      && currentTime >= setSpecificStatesCallEnableBegin-1e-13)
  {
    VLOG(1) << "-> checkStimulation";
    checkStimulation = true;

    VLOG(1) << "check if stimulation is over: duration already: " << currentTime - (lastStimulationCheckTime + 1./(setSpecificStatesCallFrequency+currentJitter))
      << ", setSpecificStatesRepeatAfterFirstCall: " << setSpecificStatesRepeatAfterFirstCall;

    // if current stimulation is over
    if (setSpecificStatesRepeatAfterFirstCall != 0
        && currentTime - (lastStimulationCheckTime + 1./(setSpecificStatesCallFrequency+currentJitter)) > setSpecificStatesRepeatAfterFirstCall)
    {
      // advance time of last call to specificStates
      LOG(DEBUG) << " old lastStimulationCheckTime: " << fiberDataCurrentPoint.lastStimulationCheckTime << ", currentJitter: " << currentJitter << ", add " << 1./(setSpecificStatesCallFrequency+currentJitter);
      fiberDataCurrentPoint.lastStimulationCheckTime += 1./(setSpecificStatesCallFrequency+currentJitter);

      LOG(DEBUG) << " new lastStimulationCheckTime: " << fiberDataCurrentPoint.lastStimulationCheckTime;

      // compute new jitter value
      double jitterFactor = 0.0;
      if (setSpecificStatesFrequencyJitter.size() > 0)
        jitterFactor = setSpecificStatesFrequencyJitter[jitterIndex % setSpecificStatesFrequencyJitter.size()];
      currentJitter = jitterFactor * setSpecificStatesCallFrequency;
      LOG(DEBUG) << " jitterIndex: " << jitterIndex << ", new jitterFactor: " << jitterFactor << ", currentJitter: " << currentJitter;
      jitterIndex++;

      checkStimulation = false;
    }
  }

  // instead of calling setSpecificStates, directly determine whether to stimulate from the firingEvents file
  int firingEventsIndex = round(currentTime * setSpecificStatesCallFrequency);

  bool stimulate =
    checkStimulation
    && firingEvents_[firingEventsIndex % firingEvents_.size()][motorUnitNo % firingEvents_[firingEventsIndex % firingEvents_.size()].size()];

  if (checkStimulation)
  {
    VLOG(1) << "setSpecificStatesCallFrequency: " << setSpecificStatesCallFrequency << ", firingEventsIndex: " << firingEventsIndex << ", fires: "
      << firingEvents_[firingEventsIndex % firingEvents_.size()][motorUnitNo % firingEvents_[firingEventsIndex % firingEvents_.size()].size()];
    VLOG(1) << "currentPointIsInCenter: " << currentPointIsInCenter;
  }

  if (stimulate)
  {
    fiberHasBeenStimulated_[fiberDataNo] = true;
  }

  if (stimulate && currentPointIsInCenter)
  {
    // if this is the first point in time of the current stimulation, log stimulation time
    if (!fiberDataCurrentPoint.currentlyStimulating)
    {
      fiberDataCurrentPoint.currentlyStimulating = true;
      Control::StimulationLogging::logStimulationBegin(currentTime, fiberDataCurrentPoint.motorUnitNo, fiberDataCurrentPoint.fiberNoGlobal);
    }

    LOG(DEBUG) << "stimulate fiber " << fiberDataCurrentPoint.fiberNoGlobal << ", MU " << motorUnitNo << " at t=" << currentTime;
    LOG(DEBUG) << "  motorUnitNo: " << motorUnitNo << " (" << motorUnitNo % firingEvents_[firingEventsIndex % firingEvents_.size()].size() << ")";
    LOG(DEBUG) << "  firing events index: " << firingEventsIndex << " (" << firingEventsIndex % firingEvents_.size() << ")";
    LOG(DEBUG) << "  setSpecificStatesCallEnableBegin: " << setSpecificStatesCallEnableBegin << ", lastStimulationCheckTime: " << lastStimulationCheckTime
      << ", stimulation already for " << + 1./(setSpecificStatesCallFrequency+currentJitter);

    LOG(INFO) << "t: " << currentTime << ", stimulate fiber " << fiberDataCurrentPoint.fiberNoGlobal << ", MU " << motorUnitNo;
  }

  if (!stimulate)
  {
    fiberDataCurrentPoint.currentlyStimulating = false;
  }

  bool stimulateCurrentPoint = stimulate && currentPointIsInCenter;
  return stimulateCurrentPoint;
}

// methods to improve speed by only computing states that are not in equilibrium
template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
equilibriumAccelerationUpdate(const Vc::double_v statesPreviousValues[], int pointBuffersNo)
{
  // check if states for the current point are at their equilibrium
  if (disableComputationWhenStatesAreCloseToEquilibrium_)
  {
    bool statesAreAtEquilibrium = true;

    // if the current point is not constant and therefore was computed
    if (fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo] != constant)
    {
      // check if it now was indeed constant
      // loop over all states
      for (int stateNo = 0; stateNo < nStates; stateNo++)
      {
        // compute relative change
        const Vc::double_v newValue = fiberPointBuffers_[pointBuffersNo].states[stateNo];
        const Vc::double_v oldValue = statesPreviousValues[stateNo];

        Vc::double_v relativeChange;
        // if any value is 0, use absolute error
        if (fabs(newValue.min()) < 1e-11)
        {
          relativeChange = Vc::abs(newValue - oldValue);
        }
        else
        {
          // values are not zero, use relative error
          relativeChange = Vc::abs((newValue - oldValue) / newValue);
        }
        double changeValue = relativeChange.max();


        // if change is higher than tolerance
        if (changeValue > 1e-5)
        {
          //LOG(INFO) << "point " << pointBuffersNo << " state " << stateNo << " " << oldValue << "->" << newValue << ": " << relativeChange << " " << changeValue;

          statesAreAtEquilibrium = false;
          break;
        }
      }
    }
/*
    // for debugging
    if (statesAreAtEquilibrium)
    {
      fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo] = constant;
    }
    else
    {
      fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo] = not_constant;
    }
*/

    // if the current point is constant but was computed because the neighbours were not constant,
    // check if this condition still holds
    if (fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo] == neighbour_not_constant
      && statesAreAtEquilibrium
    )
    {
      // check if one of the neighbours is still not constant, if it is not, set own point to constant
      if (pointBuffersNo > 0)
        if (fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo-1] != not_constant)
        {
          if (fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo] != constant)
          {
            fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo] = constant;
            nFiberPointBufferStatesCloseToEquilibrium_++;
          }
        }

      if (pointBuffersNo < fiberPointBuffers_.size()-1)
        if (fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo+1] != not_constant)
        {
          if (fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo] != constant)
          {
            fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo] = constant;
            nFiberPointBufferStatesCloseToEquilibrium_++;
          }
        }

    }

    // if own point was previously not constant, eventually set to constant
    if (fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo] == not_constant
      || fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo] == neighbour_not_constant
    )
    {
      if (statesAreAtEquilibrium)
      {
        if (fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo] != constant)
        {
          fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo] = constant;
          nFiberPointBufferStatesCloseToEquilibrium_++;
        }
      }
      else
      {
        // own point is not constant

        if (fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo] == constant)
          nFiberPointBufferStatesCloseToEquilibrium_--;

        fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo] = not_constant;

        // set neighbouring point buffers to "neighbour_not_constant"
        if (pointBuffersNo > 0)
          if (fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo-1] == constant)
          {
            if (fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo-1] == constant)
              nFiberPointBufferStatesCloseToEquilibrium_--;

            fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo-1] = neighbour_not_constant;
          }
        if (pointBuffersNo < fiberPointBuffers_.size()-1)
          if (fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo+1] == constant)
          {
            if (fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo+1] == constant)
              nFiberPointBufferStatesCloseToEquilibrium_--;
            fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo+1] = neighbour_not_constant;
          }
      }
    }
  }
}

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
bool FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
isEquilibriumAccelerationCurrentPointDisabled(bool stimulateCurrentPoint, int pointBuffersNo)
{
  if (disableComputationWhenStatesAreCloseToEquilibrium_)
  {
    if (stimulateCurrentPoint)
    {
      fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo] = not_constant;

      // set neighbouring point buffers to "neighbour_not_constant"
      if (pointBuffersNo > 0)
        if (fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo-1] == constant)
        {
          if (fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo-1] == constant)
            nFiberPointBufferStatesCloseToEquilibrium_--;
          fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo-1] = neighbour_not_constant;
        }
      if (pointBuffersNo < fiberPointBuffers_.size()-1)
        if (fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo+1] == constant)
        {
          if (fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo+1] == constant)
            nFiberPointBufferStatesCloseToEquilibrium_--;
          fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo+1] = neighbour_not_constant;
        }
    }

    // do not compute current point if the states are at the equilibrium anyways and won't change by the computation
    if (fiberPointBuffersStatesAreCloseToEquilibrium_[pointBuffersNo] == constant)
    {
      return true;
    }
  }
  return false;
}
