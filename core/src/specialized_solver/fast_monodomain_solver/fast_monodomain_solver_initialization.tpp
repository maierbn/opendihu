#include "specialized_solver/fast_monodomain_solver/fast_monodomain_solver_base.h"

#include "partition/rank_subset.h"
#include "control/diagnostic_tool/stimulation_logging.h"
#include <Vc/Vc>

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
FastMonodomainSolverBase(const DihuContext &context) :
  specificSettings_(context.getPythonConfig()), nestedSolvers_(context),
  initialized_(false)
{
  // initialize output writers
  this->outputWriterManager_.initialize(context, specificSettings_);
}

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
initialize()
{
  // only initialize once
  if (initialized_)
    return;

  if (!std::is_same<DiffusionTimeSteppingScheme,
        TimeSteppingScheme::ImplicitEuler<typename DiffusionTimeSteppingScheme::DiscretizableInTime>
      >::value
     && !std::is_same<DiffusionTimeSteppingScheme,
        TimeSteppingScheme::CrankNicolson<typename DiffusionTimeSteppingScheme::DiscretizableInTime>
      >::value)
  {
    LOG(FATAL) << "Timestepping scheme of diffusion in FastMonodomainSolver must be either ImplicitEuler or CrankNicolson!";
  }

  // add this solver to the solvers diagram, which is an ASCII art representation that will be created at the end of the simulation.
  DihuContext::solverStructureVisualizer()->addSolver("FastMonodomainSolver", true);   // hasInternalConnectionToFirstNestedSolver=true (the last argument) means slot connector data is shared with the first subsolver

  // indicate in solverStructureVisualizer that now a child solver will be initialized
  DihuContext::solverStructureVisualizer()->beginChild();

  // initialize all the nested solvers
  nestedSolvers_.initialize();

  // indicate in solverStructureVisualizer that the child solver initialization is done
  DihuContext::solverStructureVisualizer()->endChild();

  // initialize motor unit numbers and firing times
  fiberDistributionFilename_ = specificSettings_.getOptionString("fiberDistributionFile", "");
  firingTimesFilename_ = specificSettings_.getOptionString("firingTimesFile", "");
  onlyComputeIfHasBeenStimulated_ = specificSettings_.getOptionBool("onlyComputeIfHasBeenStimulated", true);
  disableComputationWhenStatesAreCloseToEquilibrium_ = specificSettings_.getOptionBool("disableComputationWhenStatesAreCloseToEquilibrium", true);
  valueForStimulatedPoint_ = specificSettings_.getOptionDouble("valueForStimulatedPoint", 20.0);

  // output warning if there are output writers
  if (this->outputWriterManager_.hasOutputWriters())
  {
    LOG(WARNING) << "The FastMonodomainSolver has output writers. Those will output all states of the 0D problem every splitting step. "
      << "This may be slow! Maybe you do not want to have this because the 1D fiber output writers can also be called.";
  }

  initializeCellMLSourceFile();
  std::shared_ptr<Partition::RankSubset> rankSubset = nestedSolvers_.data().functionSpace()->meshPartition()->rankSubset();

  LOG(DEBUG) << "config: " << specificSettings_;
  LOG(DEBUG) << "fiberDistributionFilename: " << fiberDistributionFilename_;
  LOG(DEBUG) << "firingTimesFilename: " << firingTimesFilename_;

  // parse firingTimesFilename_
  std::string firingTimesFileContents = MPIUtility::loadFile(firingTimesFilename_, rankSubset->mpiCommunicator());

  // parse file contents of firing times file
  while (!firingTimesFileContents.empty())
  {
    // extract line
    std::size_t lineEndPos = firingTimesFileContents.find("\n");
    if (lineEndPos == std::string::npos)
      break;

    std::string line = firingTimesFileContents.substr(0, lineEndPos);
    firingTimesFileContents.erase(0, lineEndPos+1);

    firingEvents_.push_back(std::vector<bool>());

    // parse line
    while (!line.empty())
    {
      int entry = atoi(line.c_str());

      firingEvents_.back().push_back((bool)(entry));

      std::size_t pos = line.find_first_of("\t ");
      if (pos == std::string::npos)
        break;
      line.erase(0, pos+1);

      while (!isdigit(line[0]) && !line.empty())
        line.erase(0,1);
    }
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
        nInstancesToCompute_ += fiberFunctionSpace->nDofsGlobal();

        assert(fiberDataNo < fiberData_.size());
        assert(fiberFunctionSpace);
        assert(motorUnitNo_.size() > 0);
        
        LOG(DEBUG) << "Fiber " << fiberNoGlobal << " (i,j)=(" << i << "," << j << ") is MU " << motorUnitNo_[fiberNoGlobal % motorUnitNo_.size()];

        fiberData_.at(fiberDataNo).valuesLength = fiberFunctionSpace->nDofsGlobal();
        fiberData_.at(fiberDataNo).fiberNoGlobal = fiberNoGlobal;
        fiberData_.at(fiberDataNo).motorUnitNo = motorUnitNo_[fiberNoGlobal % motorUnitNo_.size()];
        
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
              if (firstStimulationTime == -1 || firstStimulationTimeMotorUnit < firstStimulationTime)
              {
                // store new firstStimulationTime and save motor unit no.
                firstStimulationMotorUnitNo = motorUnitNo_[fiberNoGlobal % motorUnitNo_.size()];
                firstStimulationTime = firstStimulationTimeMotorUnit;
              }
              break;
            }
          }
        }

        // increase index for fiberData_ struct
        fiberDataNo++;
      }
    }
  }

  LOG(INFO) << "Time of first stimulation: " << firstStimulationTime << ", motor unit " << firstStimulationMotorUnitNo;

  // get the states and algebraics no.s to be transferred as slot connector data
  CellmlAdapterType &cellmlAdapter = instances[0].timeStepping1().instancesLocal()[0].discretizableInTime();
  statesForTransfer_ = cellmlAdapter.statesForTransfer();
  algebraicsForTransfer_ = cellmlAdapter.algebraicsForTransfer();


  int nInstancesLocalCellml;
  int nAlgebraicsLocalCellml;
  int nParametersPerInstance;
  cellmlAdapter.getNumbers(nInstancesLocalCellml, nAlgebraicsLocalCellml, nParametersPerInstance);

  // make parameterValues() of cellmlAdapter.data() available
  cellmlAdapter.data().prepareParameterValues();
  double *parameterValues = cellmlAdapter.data().parameterValues();   //< contains nAlgebraics parameters for all instances, in struct of array ordering (p0inst0, p0inst1, p0inst2,...)

  int nVcVectors = (nInstancesToCompute_ + Vc::double_v::Size - 1) / Vc::double_v::Size;

  fiberPointBuffers_.resize(nVcVectors);
  fiberPointBuffersAlgebraicsForTransfer_.resize(nVcVectors);
  fiberPointBuffersParameters_.resize(nVcVectors);
  fiberPointBuffersStatesAreCloseToEquilibrium_.resize(nVcVectors, not_constant);
  nFiberPointBufferStatesCloseToEquilibrium_ = 0;

  for (int i = 0; i < nVcVectors; i++)
  {
    fiberPointBuffersAlgebraicsForTransfer_[i].resize(algebraicsForTransfer_.size());
    fiberPointBuffersParameters_[i].resize(nParametersPerInstance);

    for (int j=0; j<nParametersPerInstance; j++)
    {
      fiberPointBuffersParameters_[i][j] = parameterValues[j*nAlgebraicsLocalCellml];    // note, the stride in parameterValues is "nAlgebraicsLocalCellml", not "nParametersPerInstance"

      VLOG(1) << "fiberPointBuffersParameters_ buffer no " << i << ", parameter no " << j << ", value: " <<  fiberPointBuffersParameters_[i][j];
    }
  }

  // close raw array representation of parameterValues() of cellmlAdapter.data()
  cellmlAdapter.data().restoreParameterValues();

  LOG(DEBUG) << nInstancesToCompute_ << " instances to compute, " << nVcVectors
    << " Vc vectors, size of double_v: " << Vc::double_v::Size << ", "
    << statesForTransfer_.size()-1 << " additional states for transfer, "
    << algebraicsForTransfer_.size() << " algebraics for transfer";

  // initialize state values
  for (int i = 0; i < fiberPointBuffers_.size(); i++)
  {
    if (initializeStates_ != nullptr)
    {
      initializeStates_(fiberPointBuffers_[i].states);
    }
    else
    {
      initializeStates(fiberPointBuffers_[i].states);
    }
  }

  // initialize field variable names
  for (int i = 0; i < instances.size(); i++)
  {
    std::vector<TimeSteppingScheme::Heun<CellmlAdapterType>> &innerInstances
      = instances[i].timeStepping1().instancesLocal();  // TimeSteppingScheme::Heun<CellmlAdapter...

    for (int j = 0; j < innerInstances.size(); j++, fiberNo++)
    {
      std::shared_ptr<FiberFunctionSpace> fiberFunctionSpace = innerInstances[j].data().functionSpace();
      CellmlAdapterType &cellmlAdapter = innerInstances[j].discretizableInTime();

      // loop over further states to transfer
      int furtherDataIndex = 0;
      for (int stateIndex = 1; stateIndex < statesForTransfer_.size(); stateIndex++, furtherDataIndex++)
      {
        std::string name = cellmlAdapter.data().states()->componentName(statesForTransfer_[stateIndex]);

        // get field variable
        std::vector<::Data::ComponentOfFieldVariable<FiberFunctionSpace,1>> &variable1
          = instances[i].timeStepping2().instancesLocal()[j].getSlotConnectorData()->variable1;

        if (stateIndex >= variable1.size())
        {
          LOG(WARNING) << "There are " << statesForTransfer_.size() << " statesForTransfer specified in StrangSplitting, "
            << " but the diffusion solver has only " << variable1.size() << " slots to get these states. "
            << "This means that state no. " << statesForTransfer_[stateIndex] << ", \"" << name << "\" cannot be transferred." << std::endl
            << "Maybe you need to increase \"nAdditionalFieldVariables\" in the diffusion solver (\"ImplicitEuler\" or \"CrankNicolson\") "
            << "or reduce the number of entries in \"statesForTransfer\".";
        }
        else
        {
          std::shared_ptr<FieldVariable::FieldVariable<FiberFunctionSpace,1>> fieldVariableStates
            = variable1[stateIndex].values;

          fieldVariableStates->setName(name);
        }
      }

      // loop over algebraics to transfer
      for (int algebraicIndex = 0; algebraicIndex < algebraicsForTransfer_.size(); algebraicIndex++, furtherDataIndex++)
      {
        std::string name = cellmlAdapter.data().algebraics()->componentName(algebraicsForTransfer_[algebraicIndex]);

        LOG(DEBUG) << "algebraicIndex " << algebraicIndex << ", algebraic " << algebraicsForTransfer_[algebraicIndex] << ", name: \"" << name << "\".";

        assert (i < instances.size());
        assert (j < instances[i].timeStepping2().instancesLocal().size());

        // get field variable
        std::vector<::Data::ComponentOfFieldVariable<FiberFunctionSpace,1>> &variable2
          = instances[i].timeStepping2().instancesLocal()[j].getSlotConnectorData()->variable2;

        if (algebraicIndex >= variable2.size())
        {
          LOG(WARNING) << "There are " << algebraicsForTransfer_.size() << " algebraicsForTransfer specified in StrangSplitting, "
            << " but the diffusion solver has only " << variable2.size() << " slots to get these algebraics. "
            << "This means that algebraic no. " << algebraicsForTransfer_[algebraicIndex] << ", \"" << name << "\" cannot be transferred." << std::endl
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
    }
  }

  initialized_ = true;
}

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
initializeCellMLSourceFile()
{
  // parse options
  CellmlAdapterType &cellmlAdapter = nestedSolvers_.instancesLocal()[0].timeStepping1().instancesLocal()[0].discretizableInTime();
  bool approximateExponentialFunction = cellmlAdapter.approximateExponentialFunction();

  PythonConfig specificSettingsCellML = cellmlAdapter.specificSettings();
  CellmlSourceCodeGenerator &cellmlSourceCodeGenerator = cellmlAdapter.cellmlSourceCodeGenerator();

  // determine filename of library
  std::stringstream s;
  s << "lib/"+StringUtility::extractBasename(cellmlSourceCodeGenerator.sourceFilename()) << "_fast_monodomain.so";
  std::string libraryFilename = s.str();

  //std::shared_ptr<Partition::RankSubset> rankSubset = nestedSolvers_.data().functionSpace()->meshPartition()->rankSubset();
  int ownRankNoCommWorld = DihuContext::ownRankNoCommWorld();

  if (ownRankNoCommWorld == 0)
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
    cellmlSourceCodeGenerator.generateSourceFileVcFastMonodomain(sourceToCompileFilename, approximateExponentialFunction);

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
