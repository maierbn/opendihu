#include "specialized_solver/fast_monodomain_solver/fast_monodomain_solver_base.h"

#include <vc_or_std_simd.h>  // this includes <Vc/Vc> or a Vc-emulating wrapper of <experimental/simd> if available
#include "partition/rank_subset.h"
#include "control/diagnostic_tool/stimulation_logging.h"
#include <random>

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
FastMonodomainSolverBase(const DihuContext &context) :
  specificSettings_(context.getPythonConfig()), nestedSolvers_(context),
  compute0DInstance_(nullptr), computeMonodomain_(nullptr), initializeStates_(nullptr), useVc_(true), initialized_(false)
{
  // initialize output writers
  this->outputWriterManager_.initialize(context, specificSettings_);
}

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
initialize()
{
  LOG_SCOPE_FUNCTION;

  LOG(DEBUG) << "initialize FastMonodomainSolverBase";

  // only initialize once
  if (initialized_)
    return;

  // if (!std::is_same<DiffusionTimeSteppingScheme,
  //       TimeSteppingScheme::ImplicitEuler<typename DiffusionTimeSteppingScheme::DiscretizableInTime>
  //     >::value
  //    && !std::is_same<DiffusionTimeSteppingScheme,
  //       TimeSteppingScheme::CrankNicolson<typename DiffusionTimeSteppingScheme::DiscretizableInTime>
  //     >::value)
  // {
  //   LOG(FATAL) << "Timestepping scheme of diffusion in FastMonodomainSolver must be either ImplicitEuler or CrankNicolson!";
  // }

  // disable the generation of the source code and compilation of the library in the nested CellmlAdapter
  // loop over instances of the CellML adapter
  for (int i = 0; i < nestedSolvers_.instancesLocal().size(); i++)
  {
    std::vector<TimeSteppingScheme::Heun<CellmlAdapterType>> &innerInstances
      = nestedSolvers_.instancesLocal()[i].timeStepping1().instancesLocal();  // TimeSteppingScheme::Heun<CellmlAdapter...

    for (int j = 0; j < innerInstances.size(); j++)
    {
      innerInstances[j].discretizableInTime().setCreateOwnRhsRoutine(false);
    }
  }

  // add this solver to the solvers diagram, which is an ASCII art representation that will be created at the end of the simulation.
  DihuContext::solverStructureVisualizer()->addSolver("FastMonodomainSolver", true);   // hasInternalConnectionToFirstNestedSolver=true (the last argument) means slot connector data is shared with the first subsolver

  // indicate in solverStructureVisualizer that now a child solver will be initialized
  DihuContext::solverStructureVisualizer()->beginChild();

  // initialize all the nested solvers
  nestedSolvers_.initialize();

  // indicate in solverStructureVisualizer that the child solver initialization is done
  DihuContext::solverStructureVisualizer()->endChild();

  // parse various options from the Python script
  fiberDistributionFilename_ = specificSettings_.getOptionString("fiberDistributionFile", "");
  firingTimesFilename_ = specificSettings_.getOptionString("firingTimesFile", "");
  onlyComputeIfHasBeenStimulated_ = specificSettings_.getOptionBool("onlyComputeIfHasBeenStimulated", true);
  disableComputationWhenStatesAreCloseToEquilibrium_ = specificSettings_.getOptionBool("disableComputationWhenStatesAreCloseToEquilibrium", true);
  valueForStimulatedPoint_ = specificSettings_.getOptionDouble("valueForStimulatedPoint", 20.0);
  neuromuscularJunctionRelativeSize_ = specificSettings_.getOptionDouble("neuromuscularJunctionRelativeSize", 0.0);
  generateGpuSource_ = specificSettings_.getOptionBool("generateGPUSource", true);

  // output warning if there are output writers
  if (this->outputWriterManager_.hasOutputWriters())
  {
    LOG(WARNING) << "The FastMonodomainSolver has output writers. Those will output all states of the 0D problem every splitting step. "
      << "This may be slow! Maybe you do not want to have this because the 1D fiber output writers can also be called.";
  }

  // determine optimization type and if the code should use Vc or GPU code
  CellmlAdapterType &cellmlAdapter = nestedSolvers_.instancesLocal()[0].timeStepping1().instancesLocal()[0].discretizableInTime();
  optimizationType_ = cellmlAdapter.optimizationType();
  LOG(DEBUG) << "optimizationType: \"" << optimizationType_ << "\".";

  if (optimizationType_ == "gpu")
    useVc_ = false;
  else if (optimizationType_ == "simd")
    useVc_ = false;
  else if (optimizationType_ == "vc")
    useVc_ = true;
  else
  {
    LOG(ERROR) << "FastMonodomainSolver is used with invalid \"optimizationType\": \"" << optimizationType_
      << "\". Valid options are \"vc\", \"simd\" or \"gpu\". Now using \"vc\".";
    useVc_ = true;
    optimizationType_ = "vc";
  }

  std::shared_ptr<Partition::RankSubset> rankSubset = nestedSolvers_.data().functionSpace()->meshPartition()->rankSubset();

  LOG(DEBUG) << "config: " << specificSettings_;
  LOG(DEBUG) << "fiberDistributionFilename: " << fiberDistributionFilename_;
  LOG(DEBUG) << "firingTimesFilename: " << firingTimesFilename_;

  // load the firing times of the motor units from a file
  initializeFiringTimes();

  // initialize all other internal data structures, also the data buffers used for GPU computations
  initializeDataStructures();

  if (useVc_)
  {
    // if the optimization type is "vc", create, compile, link and load the according C++-Source for CPU
    initializeCellMLSourceFileVc();
  }
  else
  {
    // if the optimzation type is GPU, create, compile, link and load the according C++-sources for GPU
    initializeCellMLSourceFileGpu();
    initializeValuesOnGpu();
  }

  // initialize state values
  for (int i = 0; i < fiberPointBuffers_.size(); i++)
  {
    // if an initialization function is given, use it to initialize the state values
    if (initializeStates_ != nullptr)
    {
      initializeStates_(fiberPointBuffers_[i].states);
    }
    else
    {
      initializeStates(fiberPointBuffers_[i].states);
    }
  }
  setComputeStateInformation_ = false;

  // initialize the variable names where field variables are connector via connector slots
  initializeFieldVariableNames();

  initialized_ = true;
}

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
initializeFiringTimes()
{
  std::shared_ptr<Partition::RankSubset> rankSubset = nestedSolvers_.data().functionSpace()->meshPartition()->rankSubset();

  // parse firingTimesFilename_
  std::string firingTimesFileContents = MPIUtility::loadFile(firingTimesFilename_, rankSubset->mpiCommunicator());
  gpuFiringEventsNColumns_ = 0;
  gpuFiringEventsNRows_ = 0;

  // parse file contents of firing times file, loop over rows
  while (!firingTimesFileContents.empty())
  {
    // determine end of file
    std::size_t lineEndPos = firingTimesFileContents.find("\n");
    if (lineEndPos == std::string::npos)
      break;

    // extract line from file contents
    std::string line = firingTimesFileContents.substr(0, lineEndPos);
    firingTimesFileContents.erase(0, lineEndPos+1);

    //variable has the following layout: firingEvents_[timeStepNo][motorUnitNo]
    firingEvents_.push_back(std::vector<bool>());

    // parse line, loop over columns
    int columnNo = 0;
    for (; !line.empty(); columnNo++)
    {
      int entry = atoi(line.c_str());

      // the variable firingEvents_ is for normal execution
      firingEvents_.back().push_back((bool)(entry));

      // the variable gpuFiringEvents_ is to be send to gpu
      if (gpuFiringEventsNColumns_ == 0 || columnNo < gpuFiringEventsNColumns_)
      {
        gpuFiringEvents_.push_back(entry? 1 : 0);
      }

      // remove separator, either space or tab
      std::size_t pos = line.find_first_of("\t ");
      if (pos == std::string::npos)
      {
        if (gpuFiringEventsNColumns_ == 0)
          gpuFiringEventsNColumns_ = columnNo+1;
        break;
      }
      line.erase(0, pos+1);

      // remove all following non-digit characters
      while (!isdigit(line[0]) && !line.empty())
        line.erase(0,1);
    }
    if (gpuFiringEventsNColumns_ == 0)
    {
      gpuFiringEventsNColumns_ = columnNo;
      LOG(DEBUG) << "firing events file contains " << gpuFiringEventsNColumns_ << " columns.";
    }
    gpuFiringEventsNRows_++;
  }

  // parse fiberDistributionFile
  std::string fiberDistributionFileContents = MPIUtility::loadFile(fiberDistributionFilename_, rankSubset->mpiCommunicator());

  // parse file contents
  while (!fiberDistributionFileContents.empty())
  {
    int motorUnitNo = atoi(fiberDistributionFileContents.c_str());
    motorUnitNo_.push_back(motorUnitNo);

    std::size_t pos = fiberDistributionFileContents.find_first_of("\t ");

    if (pos == std::string::npos)
      break;

    fiberDistributionFileContents.erase(0, pos+1);
  }

  LOG(DEBUG) << "firingEvents.size: " << firingEvents_.size();
  LOG(DEBUG) << "firingEvents_:" << firingEvents_;
  LOG(DEBUG) << "motorUnitNo_: " << motorUnitNo_;

  if (motorUnitNo_.empty())
    LOG(FATAL) << "Could not parse motor units.";

  if (firingEvents_.empty())
    LOG(FATAL) << "Could not parse firing times.";
}

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
initializeDataStructures()
{
  // initialize data structures
  std::vector<typename NestedSolversType::TimeSteppingSchemeType> &instances = nestedSolvers_.instancesLocal();

  // determine number of fibers to compute on the current rank
  int nFibers = 0;
  int fiberNo = 0;
  nFibersToCompute_ = 0;

  LOG(DEBUG) << "initialize " << instances.size() << " outer instances";

  for (int i = 0; i < instances.size(); i++)
  {
    std::vector<TimeSteppingScheme::Heun<CellmlAdapterType>> &innerInstances
      = instances[i].timeStepping1().instancesLocal();  // TimeSteppingScheme::Heun<CellmlAdapter...
    nFibers += innerInstances.size();

    for (int j = 0; j < innerInstances.size(); j++, fiberNo++)
    {
      std::shared_ptr<FiberFunctionSpace> fiberFunctionSpace = innerInstances[j].data().functionSpace();
      std::shared_ptr<Partition::RankSubset> rankSubset = fiberFunctionSpace->meshPartition()->rankSubset();
      int computingRank = fiberNo % rankSubset->size();

      LOG(DEBUG) << "instance (inner,outer)=(i,j)=(" << i << "," << j << ")/(" << instances.size() << "," << innerInstances.size() << ")"
        << ", fiberNo " << fiberNo << ", rankSubset: " << *rankSubset << ", mesh" << fiberFunctionSpace->meshName() << ", computingRank " << computingRank << ", own rank: " << rankSubset->ownRankNo() << "/" << rankSubset->size();

      if (computingRank == rankSubset->ownRankNo())
      {
        LOG(DEBUG) << "compute (i,j)=(" << i << "," << j << "), computingRank " << computingRank
          << ", fiberNo: " << fiberNo << ", rankSubset->size(): " << rankSubset->size() << ", own: " << rankSubset->ownRankNo();
        nFibersToCompute_++;
      }
    }
  }

  fiberData_.resize(nFibersToCompute_);
  fiberHasBeenStimulated_.resize(nFibersToCompute_, false);
  LOG(DEBUG) << "nFibers: " << nFibers << ", nFibersToCompute_: " << nFibersToCompute_;

  // initialize random generator that is used to determine the stimulation point of the fiber
  std::random_device randomDevice;              // Will be used to obtain a seed for the random number engine
  std::mt19937 randomGenerator(randomDevice()); // Standard mersenne_twister_engine seeded with randomDevice()
  std::uniform_real_distribution<> randomDistribution(0.5-neuromuscularJunctionRelativeSize_/2., 0.5+neuromuscularJunctionRelativeSize_/2.);
  // now, a call to randomDistribution(randomGenerator)) will produce a uniformly distributed value between 0.5-a and 0.5+a. 0.5 corresponds to the center of the fiber.

  // determine total number of CellML instances to compute on this rank
  double firstStimulationTime = -1;
  int firstStimulationMotorUnitNo = 0;
  nInstancesToCompute_ = 0;
  int fiberDataNo = 0;
  fiberNo = 0;
  for (int i = 0; i < instances.size(); i++)
  {
    std::vector<TimeSteppingScheme::Heun<CellmlAdapterType>> &innerInstances
      = instances[i].timeStepping1().instancesLocal();  // TimeSteppingScheme::Heun<CellmlAdapter...

    for (int j = 0; j < innerInstances.size(); j++, fiberNo++)
    {
      CellmlAdapterType &cellmlAdapter = innerInstances[j].discretizableInTime();
      int fiberNoGlobal = PythonUtility::convertFromPython<int>::get(cellmlAdapter.pySetFunctionAdditionalParameter_);

      // initialize filename in stimulation logging class from current settings
      Control::StimulationLogging logging(cellmlAdapter.specificSettings_);

      std::shared_ptr<FiberFunctionSpace> fiberFunctionSpace = innerInstances[j].data().functionSpace();

      std::shared_ptr<Partition::RankSubset> rankSubset = fiberFunctionSpace->meshPartition()->rankSubset();
      int computingRank = fiberNo % rankSubset->size();

      if (computingRank == rankSubset->ownRankNo())
      {
        LOG(DEBUG) << "compute (i,j)=(" << i << "," << j << "), computingRank " << computingRank
          << ", fiberNo: " << fiberNo << ", fiberDataNo: " << fiberDataNo
          << ", rankSubset->size(): " << rankSubset->size() << ", own: " << rankSubset->ownRankNo();
        nInstancesToComputePerFiber_ = fiberFunctionSpace->nDofsGlobal();
        nInstancesToCompute_ += nInstancesToComputePerFiber_;

        assert(fiberDataNo < fiberData_.size());
        assert(fiberFunctionSpace);
        assert(motorUnitNo_.size() > 0);

        LOG(DEBUG) << "Fiber " << fiberNoGlobal << " (i,j)=(" << i << "," << j << ") is MU " << motorUnitNo_[fiberNoGlobal % motorUnitNo_.size()];

        fiberData_.at(fiberDataNo).valuesLength = nInstancesToComputePerFiber_;
        fiberData_.at(fiberDataNo).fiberNoGlobal = fiberNoGlobal;
        fiberData_.at(fiberDataNo).motorUnitNo = motorUnitNo_[fiberNoGlobal % motorUnitNo_.size()];

        // determine neuromuscular junction position, it is offset by a random value from the center
        if (neuromuscularJunctionRelativeSize_ <= 0.0)
        {
          fiberData_.at(fiberDataNo).fiberStimulationPointIndex = (int)(fiberData_.at(fiberDataNo).valuesLength / 2);
        }
        else
        {
          fiberData_.at(fiberDataNo).fiberStimulationPointIndex = (int)(fiberData_.at(fiberDataNo).valuesLength * randomDistribution(randomGenerator));
        }

        // copy settings
        fiberData_.at(fiberDataNo).setSpecificStatesCallFrequency = cellmlAdapter.setSpecificStatesCallFrequency_;
        fiberData_.at(fiberDataNo).setSpecificStatesFrequencyJitter = cellmlAdapter.setSpecificStatesFrequencyJitter_;
        fiberData_.at(fiberDataNo).setSpecificStatesRepeatAfterFirstCall = cellmlAdapter.setSpecificStatesRepeatAfterFirstCall_;
        fiberData_.at(fiberDataNo).setSpecificStatesCallEnableBegin = cellmlAdapter.setSpecificStatesCallEnableBegin_;

        fiberData_.at(fiberDataNo).currentJitter = 0;
        fiberData_.at(fiberDataNo).jitterIndex = 0;
        fiberData_.at(fiberDataNo).lastStimulationCheckTime = fiberData_.at(fiberDataNo).setSpecificStatesCallEnableBegin - 1e-13 - 1./(fiberData_.at(fiberDataNo).setSpecificStatesCallFrequency+fiberData_.at(fiberDataNo).currentJitter);

        fiberData_.at(fiberDataNo).valuesOffset = 0;
        fiberData_.at(fiberDataNo).currentlyStimulating = false;
        if (fiberDataNo > 0)
        {
          fiberData_.at(fiberDataNo).valuesOffset = fiberData_.at(fiberDataNo-1).valuesOffset + fiberData_.at(fiberDataNo-1).valuesLength;
        }

        // find out first stimulation time of any fiber
        int firingEventsIndex = round(fiberData_.at(fiberDataNo).setSpecificStatesCallEnableBegin * fiberData_.at(fiberDataNo).setSpecificStatesCallFrequency);
        if (firingEventsIndex < 0)
          firingEventsIndex = 0;

        // only if there is a chance that the current fiber will stimulate before the currently firstStimulationTime, because of setSpecificStatesCallEnableBegin
        if (firstStimulationTime == -1 || firstStimulationTime > fiberData_.at(fiberDataNo).setSpecificStatesCallEnableBegin)
        {
          // start at index in firingEvents_, that comes as setSpecificStatesCallEnableBegin, then step over next timesteps until the next stimulation is found
          int nFiringEvents = firingEvents_.size();
          for (int i = firingEventsIndex % nFiringEvents; i < nFiringEvents; i++)
          {
            // if there is a stimulation at the current timestep
            int motorUnitNo = motorUnitNo_[fiberNoGlobal % motorUnitNo_.size()];
            if (firingEvents_[i][motorUnitNo % firingEvents_[i].size()])
            {
              // compute current time
              int firstFiringEventIndex = int(firingEventsIndex / nFiringEvents)*nFiringEvents + i;
              double firstStimulationTimeMotorUnit = firstFiringEventIndex / fiberData_.at(fiberDataNo).setSpecificStatesCallFrequency;

              LOG(DEBUG) << "  Motor unit " << motorUnitNo << " fires at " << firstStimulationTimeMotorUnit;  // this output does not get called for all motor units!

              // if this time is smaller than currently saved firstStimulationTime, or firstStimulationTime has not yet been initialized
              if (fiberData_.at(fiberDataNo).setSpecificStatesCallFrequency > 1e-12 &&
                  (firstStimulationTime == -1 || firstStimulationTimeMotorUnit < firstStimulationTime))
              {
                // store new firstStimulationTime and save motor unit no.
                firstStimulationMotorUnitNo = motorUnitNo_[fiberNoGlobal % motorUnitNo_.size()];
                firstStimulationTime = firstStimulationTimeMotorUnit;
              }
              break;
            }
          }
        }

        // add values to gpuSetSpecificStatesFrequencyJitter_
        gpuSetSpecificStatesFrequencyJitter_.insert(gpuSetSpecificStatesFrequencyJitter_.end(),
                                                    cellmlAdapter.setSpecificStatesFrequencyJitter_.begin(),
                                                    cellmlAdapter.setSpecificStatesFrequencyJitter_.end());
        gpuFrequencyJitterNColumns_ = cellmlAdapter.setSpecificStatesFrequencyJitter_.size();

        // increase index for fiberData_ struct
        fiberDataNo++;
      }
    }
  }

  if (firstStimulationTime == -1)
  {
    LOG(INFO) << "No stimulation by setSpecificStates";
  }
  else
  {
    LOG(INFO) << "Time of first stimulation: " << firstStimulationTime << ", motor unit " << firstStimulationMotorUnitNo;
  }

  CellmlAdapterType &cellmlAdapter = nestedSolvers_.instancesLocal()[0].timeStepping1().instancesLocal()[0].discretizableInTime();

  // get the states and algebraics no.s to be transferred as slot connector data
  statesForTransferIndices_ = cellmlAdapter.statesForTransfer();
  algebraicsForTransferIndices_ = cellmlAdapter.algebraicsForTransfer();

  int nInstancesLocalCellml;
  int nAlgebraicsLocalCellml;
  cellmlAdapter.getNumbers(nInstancesLocalCellml, nAlgebraicsLocalCellml, nParametersPerInstance_);

  // make parameterValues() of cellmlAdapter.data() available
  cellmlAdapter.data().prepareParameterValues();
  double *parameterValues = cellmlAdapter.data().parameterValues();   //< contains nAlgebraics parameters for all instances, in struct of array ordering (p0inst0, p0inst1, p0inst2,...)

  int nVcVectors = (nInstancesToCompute_ + Vc::double_v::size() - 1) / Vc::double_v::size();

  if (useVc_)
  {
    fiberPointBuffers_.resize(nVcVectors);
    fiberPointBuffersAlgebraicsForTransfer_.resize(nVcVectors);
    fiberPointBuffersParameters_.resize(nVcVectors);
    fiberPointBuffersStatesAreCloseToEquilibrium_.resize(nVcVectors, active);
    nFiberPointBufferStatesCloseToEquilibrium_ = 0;

    for (int i = 0; i < nVcVectors; i++)
    {
      fiberPointBuffersAlgebraicsForTransfer_[i].resize(algebraicsForTransferIndices_.size());
      fiberPointBuffersParameters_[i].resize(nParametersPerInstance_);

      for (int parameterNo = 0; parameterNo < nParametersPerInstance_; parameterNo++)
      {
        fiberPointBuffersParameters_[i][parameterNo] = parameterValues[parameterNo*nAlgebraicsLocalCellml];    // note, the stride in parameterValues is "nAlgebraicsLocalCellml", not "nParametersPerInstance_"

#ifndef HAVE_STDSIMD
        VLOG(1) << "fiberPointBuffersParameters_ buffer no " << i << ", parameter no " << parameterNo
          << ", value: " <<  fiberPointBuffersParameters_[i][parameterNo];
#endif
      }
    }
  }
  else
  {
    gpuParameters_.resize(nInstancesToCompute_*nParametersPerInstance_);
    gpuAlgebraicsForTransfer_.resize(nInstancesToCompute_*algebraicsForTransferIndices_.size());
    gpuStatesForTransfer_.resize(nInstancesToCompute_*statesForTransferIndices_.size());
    gpuFiberIsCurrentlyStimulated_.resize(nFibersToCompute_, 0);

    int nElementsOnFiber = nInstancesToComputePerFiber_-1;
    gpuElementLengths_.resize(nElementsOnFiber*nFibersToCompute_);

    for (int parameterNo = 0; parameterNo < nParametersPerInstance_; parameterNo++)
    {
      for (int instanceNoToCompute = 0; instanceNoToCompute < nInstancesToCompute_; instanceNoToCompute++)
      {
        gpuParameters_[parameterNo*nInstancesToCompute_ + instanceNoToCompute] = parameterValues[parameterNo*nAlgebraicsLocalCellml];    // note, the stride in parameterValues is "nAlgebraicsLocalCellml", not "nParametersPerInstance_"
      }
    }

    gpuVmValues_.resize(nInstancesToCompute_);
    gpuMotorUnitNo_.resize(nFibersToCompute_);
    gpuFiberStimulationPointIndex_.resize(nFibersToCompute_);
    gpuLastStimulationCheckTime_.resize(nFibersToCompute_);
    gpuSetSpecificStatesCallFrequency_.resize(nFibersToCompute_);
    gpuSetSpecificStatesRepeatAfterFirstCall_.resize(nFibersToCompute_);
    gpuSetSpecificStatesCallEnableBegin_.resize(nFibersToCompute_);
    gpuCurrentJitter_.resize(nFibersToCompute_, 0.0);
    gpuJitterIndex_.resize(nFibersToCompute_, 0);

    for (int fiberDataNo = 0; fiberDataNo < nFibersToCompute_; fiberDataNo++)
    {
      gpuMotorUnitNo_[fiberDataNo]                           = fiberData_.at(fiberDataNo).motorUnitNo;
      gpuFiberStimulationPointIndex_[fiberDataNo]            = fiberData_.at(fiberDataNo).fiberStimulationPointIndex;
      gpuLastStimulationCheckTime_[fiberDataNo]              = fiberData_.at(fiberDataNo).lastStimulationCheckTime;
      gpuSetSpecificStatesCallFrequency_[fiberDataNo]        = fiberData_.at(fiberDataNo).setSpecificStatesCallFrequency;
      gpuSetSpecificStatesRepeatAfterFirstCall_[fiberDataNo] = fiberData_.at(fiberDataNo).setSpecificStatesRepeatAfterFirstCall;
      gpuSetSpecificStatesCallEnableBegin_[fiberDataNo]      = fiberData_.at(fiberDataNo).setSpecificStatesCallEnableBegin;
    }
  }

  // close raw array representation of parameterValues() of cellmlAdapter.data()
  cellmlAdapter.data().restoreParameterValues();

  LOG(DEBUG) << nInstancesToCompute_ << " instances to compute, " << nVcVectors
    << " Vc vectors, size of double_v: " << Vc::double_v::size() << ", "
    << statesForTransferIndices_.size()-1 << " additional states for transfer, "
    << algebraicsForTransferIndices_.size() << " algebraics for transfer";
}

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
initializeFieldVariableNames()
{
  // initialize data structures
  std::vector<typename NestedSolversType::TimeSteppingSchemeType> &instances = nestedSolvers_.instancesLocal();

  // initialize field variable names
  for (int i = 0; i < instances.size(); i++)
  {
    std::vector<TimeSteppingScheme::Heun<CellmlAdapterType>> &innerInstances
      = instances[i].timeStepping1().instancesLocal();  // TimeSteppingScheme::Heun<CellmlAdapter...

    for (int j = 0; j < innerInstances.size(); j++)
    {
      std::shared_ptr<FiberFunctionSpace> fiberFunctionSpace = innerInstances[j].data().functionSpace();
      CellmlAdapterType &cellmlAdapter = innerInstances[j].discretizableInTime();

      // loop over further states to transfer
      int furtherDataIndex = 0;
      for (int stateIndex = 1; stateIndex < statesForTransferIndices_.size(); stateIndex++, furtherDataIndex++)
      {
        std::string name = cellmlAdapter.data().states()->componentName(statesForTransferIndices_[stateIndex]);

        // get field variable
        std::vector<::Data::ComponentOfFieldVariable<FiberFunctionSpace,1>> &variable1
          = instances[i].timeStepping2().instancesLocal()[j].getSlotConnectorData()->variable1;

        if (stateIndex >= variable1.size())
        {
          LOG(WARNING) << "There are " << statesForTransferIndices_.size() << " statesForTransfer specified in StrangSplitting, "
            << " but the diffusion solver has only " << variable1.size() << " slots in the first variable to get these states. "
            << "This means that state no. " << statesForTransferIndices_[stateIndex] << ", \"" << name << "\" cannot be transferred." << std::endl
            << "Maybe you need to reduce the number of entries in \"statesForTransfer\". "
            << "(Note that the FastMonodomainSolver can only map algebraics to additionalFieldVariables, not states.)";
        }
        else
        {
          std::shared_ptr<FieldVariable::FieldVariable<FiberFunctionSpace,1>> fieldVariableStates
            = variable1[stateIndex].values;

          fieldVariableStates->setName(name);
        }
      }

      // get field variable
      std::vector<::Data::ComponentOfFieldVariable<FiberFunctionSpace,1>> &variable2
        = instances[i].timeStepping2().instancesLocal()[j].getSlotConnectorData()->variable2;

      // loop over algebraics to transfer
      for (int algebraicIndex = 0; algebraicIndex < algebraicsForTransferIndices_.size(); algebraicIndex++, furtherDataIndex++)
      {
        std::string name = cellmlAdapter.data().algebraics()->componentName(algebraicsForTransferIndices_[algebraicIndex]);

        LOG(DEBUG) << "algebraicIndex " << algebraicIndex << ", algebraic " << algebraicsForTransferIndices_[algebraicIndex] << ", name: \"" << name << "\".";

        assert (i < instances.size());
        assert (j < instances[i].timeStepping2().instancesLocal().size());

        if (algebraicIndex >= variable2.size())
        {
          LOG(WARNING) << "There are " << algebraicsForTransferIndices_.size() << " algebraicsForTransfer specified in StrangSplitting, "
            << " but the diffusion solver has only " << variable2.size() << " slots to get these algebraics. "
            << "This means that algebraic no. " << algebraicsForTransferIndices_[algebraicIndex] << ", \"" << name << "\" cannot be transferred." << std::endl
            << "Maybe you need to increase \"nAdditionalFieldVariables\" in  the diffusion solver (\"ImplicitEuler\" or \"CrankNicolson\") "
            << "or reduce the number of entries in \"algebraicsForTransfer\".";
        }
        else
        {
          std::shared_ptr<FieldVariable::FieldVariable<FiberFunctionSpace,1>> fieldVariableAlgebraics
            = variable2[algebraicIndex].values;

          fieldVariableAlgebraics->setName(name);
        }
      }

      // prepare extra slot for compute state information (i.e., inactive, neighbour_is_active or active)
      if (variable2.size() > algebraicsForTransferIndices_.size())
      {
        int algebraicIndex = algebraicsForTransferIndices_.size();

        // get algebraics field variable
        std::vector<::Data::ComponentOfFieldVariable<FiberFunctionSpace,1>> &variable2
          = instances[i].timeStepping2().instancesLocal()[j].getSlotConnectorData()->variable2;
        std::shared_ptr<FieldVariable::FieldVariable<FiberFunctionSpace,1>> fieldVariableAlgebraics
            = variable2[algebraicIndex].values;

        std::string name = "computeStateInformation";
        fieldVariableAlgebraics->setName(name);
        setComputeStateInformation_ = true;
      }
    }
  }
}


template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
initializeCellMLSourceFileVc()
{
  // parse options
  CellmlAdapterType &cellmlAdapter = nestedSolvers_.instancesLocal()[0].timeStepping1().instancesLocal()[0].discretizableInTime();
  bool approximateExponentialFunction = cellmlAdapter.approximateExponentialFunction();

  PythonConfig specificSettingsCellML = cellmlAdapter.specificSettings();
  CellmlSourceCodeGenerator &cellmlSourceCodeGenerator = cellmlAdapter.cellmlSourceCodeGenerator();

  std::string libraryFilename;

  // if option "libraryFilename" is given, do not create new C++ source code, use the existing library instead
  if (specificSettingsCellML.hasKey("libraryFilename"))
  {
    libraryFilename = specificSettingsCellML.getOptionString("libraryFilename", "lib.so");

    LOG(INFO) << "Loading existing library \"" << libraryFilename << "\".";
  }
  else
  {
    // option "libraryFilename" was not given, create source code for GPU and and compile it to the shared library

    // determine filename of library
    std::stringstream s;
    s << "lib/"+StringUtility::extractBasename(cellmlSourceCodeGenerator.sourceFilename()) << "_fast_monodomain.so";
    libraryFilename = s.str();

    //std::shared_ptr<Partition::RankSubset> rankSubset = nestedSolvers_.data().functionSpace()->meshPartition()->rankSubset();
    int ownRankNo = DihuContext::partitionManager()->rankSubsetForCollectiveOperations()->ownRankNo();

    if (ownRankNo == 0)
    {
      // initialize generated source code of cellml model

      // compile source file to a library

      s.str("");
      s << "src/"+StringUtility::extractBasename(cellmlSourceCodeGenerator.sourceFilename()) << "_fast_monodomain"
        << cellmlSourceCodeGenerator.sourceFileSuffix();
      std::string sourceToCompileFilename = s.str();

      // create path of library filename if it does not exist
      if (libraryFilename.find("/") != std::string::npos)
      {
        std::string path = libraryFilename.substr(0, libraryFilename.rfind("/"));
        // if directory does not yet exist, create it
        struct stat info;
        if (stat(path.c_str(), &info) != 0)
        {
          int ret = system((std::string("mkdir -p ")+path).c_str());

          if (ret != 0)
          {
            LOG(ERROR) << "Could not create path \"" << path << "\".";
          }
        }
      }

      // generate library

      LOG(DEBUG) << "generate source file \"" << sourceToCompileFilename << "\".";

      // create source file
      cellmlSourceCodeGenerator.generateSourceFileFastMonodomain(sourceToCompileFilename, approximateExponentialFunction);

      // create path for library file
      if (libraryFilename.find("/") != std::string::npos)
      {
        std::string path = libraryFilename.substr(0, libraryFilename.rfind("/"));
        int ret = system((std::string("mkdir -p ")+path).c_str());

        if (ret != 0)
        {
          LOG(ERROR) << "Could not create path \"" << path << "\" for library file.";
        }
      }

      std::stringstream compileCommand;

      // load compiler flags
      std::string compilerFlags = specificSettingsCellML.getOptionString("compilerFlags", "-O3 -march=native -fPIC -finstrument-functions -ftree-vectorize -fopt-info-vec-optimized=vectorizer_optimized.log -shared ");

    #ifdef NDEBUG
      if (compilerFlags.find("-O3") == std::string::npos)
      {
        LOG(WARNING) << "\"compilerFlags\" does not contain \"-O3\", this may be slow.";
      }
      if (compilerFlags.find("-m") == std::string::npos)
      {
        LOG(WARNING) << "\"compilerFlags\" does not contain any \"-m\" flag, such as \"-march=native\". "
          << "Make sure that SIMD instructions sets (SEE, AVX-2 etc.) are the same in opendihu and the compiled library. \n"
          << " If unsure, use \"-O3 -march-native\".";
      }
    #endif
      // for GPU: -ta=host,tesla,time

      // compose compile command
      std::stringstream s;
      s << cellmlSourceCodeGenerator.compilerCommand() << " " << sourceToCompileFilename << " "
        << compilerFlags << " " << cellmlSourceCodeGenerator.additionalCompileFlags() << " ";

      std::string compileCommandOptions = s.str();

      compileCommand << compileCommandOptions
        << " -o " << libraryFilename;

      int ret = system(compileCommand.str().c_str());
      if (ret != 0)
      {
        LOG(ERROR) << "Compilation failed. Command: \"" << compileCommand.str() << "\".";

        // remove "-fopenmp" in the compile command
        std::string newCompileCommand = compileCommand.str();
        std::string strToReplace = "-fopenmp";
        std::size_t pos = newCompileCommand.find(strToReplace);
        newCompileCommand.replace(pos, strToReplace.length(), "");

        // remove -foffload="..."strToReplace = "-fopenmp";
        pos = newCompileCommand.find("-foffload=\"");
        std::size_t pos2 = newCompileCommand.find("\"", pos+11);
        newCompileCommand.replace(pos, pos2-pos+1, "");

        LOG(INFO) << "Retry without offloading, command: \n" << newCompileCommand;

        // execute new compilation command
        int ret = system(newCompileCommand.c_str());
        if (ret != 0)
        {
          LOG(ERROR) << "Compilation failed again.";
        }
        else
        {
          LOG(DEBUG) << "Compilation successful.";
        }
      }
      else
      {
        LOG(DEBUG) << "Compilation successful. Command: \"" << compileCommand.str() << "\".";
      }
    }

    // barrier disabled because in interferes with the barrier in 00_source_code_generator_base.cpp
    //LOG(ERROR) << "MPI barrier in fast_monodomain_solver on MPI_COMM_WORLD";

    // wait on all ranks until conversion is finished
    MPIUtility::handleReturnValue(MPI_Barrier(DihuContext::partitionManager()->rankSubsetForCollectiveOperations()->mpiCommunicator()), "MPI_Barrier");
  }

  // load the rhs library
  void *handle = CellmlAdapterType::loadRhsLibraryGetHandle(libraryFilename);

  compute0DInstance_ = (void (*)(Vc::double_v [], std::vector<Vc::double_v> &, double, double, bool, bool, std::vector<Vc::double_v> &, const std::vector<int> &, double)) dlsym(handle, "compute0DInstance");
  initializeStates_ = (void (*)(Vc::double_v states[])) dlsym(handle, "initializeStates");

  LOG(DEBUG) << "compute0DInstance_: " << (compute0DInstance_==nullptr? "no" : "yes") << ", initializeStates_: " << (initializeStates_==nullptr? "no" : "yes");

  if (compute0DInstance_ == nullptr || initializeStates_ == nullptr)
  {
    LOG(FATAL) << "Could not load functions from library \"" << libraryFilename << "\".";
  }
}
