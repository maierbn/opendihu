#include "specialized_solver/fast_monodomain_solver/fast_monodomain_solver_base.h"

#include "partition/rank_subset.h"
#include "control/diagnostic_tool/stimulation_logging.h"
#include <Vc/Vc>
#include <random>

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
initializeCellMLSourceFileGpu()
{
  // parse options
  CellmlAdapterType &cellmlAdapter = nestedSolvers_.instancesLocal()[0].timeStepping1().instancesLocal()[0].discretizableInTime();
  bool approximateExponentialFunction = cellmlAdapter.approximateExponentialFunction();

  PythonConfig specificSettingsCellML = cellmlAdapter.specificSettings();
  CellmlSourceCodeGenerator &cellmlSourceCodeGenerator = cellmlAdapter.cellmlSourceCodeGenerator();

  // create source code for the rhs part
  std::string headerCode;
  std::string mainCode;
  const bool hasAlgebraicsForTransfer = !algebraicsForTransferIndices_.empty();
  cellmlSourceCodeGenerator.generateSourceFastMonodomainGpu(approximateExponentialFunction,
                                                            1, nInstancesToComputePerFiber_, nParametersPerInstance_,
                                                            hasAlgebraicsForTransfer,
                                                            headerCode, mainCode);

  int ownRankNoCommWorld = DihuContext::ownRankNoCommWorld();

  // determine filename of library
  std::stringstream s;
  s << "lib/" << StringUtility::extractBasename(cellmlSourceCodeGenerator.sourceFilename())
    << "_gpu_fast_monodomain." << ownRankNoCommWorld << ".so";
  std::string libraryFilename = s.str();

  s.str("");
  s << "src/" << StringUtility::extractBasename(cellmlSourceCodeGenerator.sourceFilename())
    << "_gpu_fast_monodomain." << ownRankNoCommWorld << ".cpp";
  std::string sourceToCompileFilename = s.str();

  // generate library
  LOG(DEBUG) << "initializeCellMLSourceFileGpu: generate source file \"" << sourceToCompileFilename << "\".";

  if (generateGpuSource_)
  {
    generateMonodomainSolverGpuSource(sourceToCompileFilename, headerCode, mainCode);
  }
  else
  {
    LOG(WARNING) << "In FastMonodomainSolver for GPU: \"generateGpuSource\" is set to False, i.e. no code will be generated."
      << " Instead, the existing source \"" << sourceToCompileFilename << "\" will be compiled.";
  }

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

  // compose compile command
  s.str("");
  s << cellmlSourceCodeGenerator.compilerCommand() << " " << sourceToCompileFilename << " "
    << compilerFlags << " " << cellmlSourceCodeGenerator.additionalCompileFlags() << " ";

  std::string compileCommandOptions = s.str();

  std::stringstream compileCommand;
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

  // wait on all ranks until compilation is finished
  MPIUtility::handleReturnValue(MPI_Barrier(DihuContext::partitionManager()->rankSubsetForCollectiveOperations()->mpiCommunicator()), "MPI_Barrier");

  // load the rhs library
  void *handle = CellmlAdapterType::loadRhsLibraryGetHandle(libraryFilename);

  computeMonodomain_ = (void (*)(const double *parameters,
                                double *algebraicsForTransfer, double *statesForTransfer, const double *elementLengths,
                                double startTime, double timeStepWidthSplitting, int nTimeStepsSplitting, double dt0D, int nTimeSteps0D, double dt1D, int nTimeSteps1D,
                                double prefactor, double valueForStimulatedPoint)) dlsym(handle, "computeMonodomain");

  initializeArrays_ = (void (*)(const double *statesOneInstance, const int *algebraicsForTransferIndicesParameter, const int *statesForTransferIndicesParameter,
                                const char *firingEventsParameter, const double *setSpecificStatesFrequencyJitterParameter, const int *motorUnitNoParameter,
                                const int *fiberStimulationPointIndexParameter, const double *lastStimulationCheckTimeParameter,
                                const double *setSpecificStatesCallFrequencyParameter, const double *setSpecificStatesRepeatAfterFirstCallParameter,
                                const double *setSpecificStatesCallEnableBeginParameter)) dlsym(handle, "initializeArrays");

  LOG(DEBUG) << "computeMonodomain: " << (computeMonodomain_==nullptr? "no" : "yes")
    << ", initializeArrays: " << (initializeArrays_==nullptr? "no" : "yes");

  if (computeMonodomain_ == nullptr || initializeArrays_ == nullptr)
  {
    LOG(INFO) << "Compile command:\n" << compileCommand.str();
    LOG(FATAL) << "Could not load functions from library \"" << libraryFilename << "\".";
  }
}

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
generateMonodomainSolverGpuSource(std::string outputFilename, std::string headerCode, std::string mainCode)
{
  LOG(TRACE) << "addMonodomainSolverGpuSource";

  // here, filename contains already c++ code generated from the cellML model that computes the 0D terms
  // we now add code to solve the whole monodomain equation

  std::stringstream sourceCode;
  sourceCode << headerCode;
  
  int nFibersToCompute = 1;  // nFibersToCompute_
  int nInstancesToCompute = nFibersToCompute*nInstancesToComputePerFiber_;

  sourceCode << R"(
#pragma omp end declare target

// global size constants
const int nInstancesPerFiber = )" << nInstancesToComputePerFiber_ << R"(;
const int nElementsOnFiber = )" << nInstancesToComputePerFiber_-1 << R"(;
const int nFibersToCompute = )" << nFibersToCompute << R"(;
const long long nInstancesToCompute = )" << nInstancesToCompute << R"(;  // = nInstancesPerFiber*nFibersToCompute;
const int nStates = )" << nStates << R"(;
const int nAlgebraics = )" << nAlgebraics << R"(;
const int firingEventsNRows = )" << gpuFiringEventsNRows_ << R"(;
const int firingEventsNColumns = )" << gpuFiringEventsNColumns_ << R"(;
const int frequencyJitterNColumns = )" << gpuFrequencyJitterNColumns_ << R"(;
const int nStatesTotal = )" << nInstancesToCompute*nStates << R"(;  // = nInstancesToCompute*nStates;
const int nAlgebraicsTotal = )" << nInstancesToCompute*nAlgebraics << R"(;  // = nInstancesToCompute*nAlgebraics;
const int nParametersTotal = )" << nInstancesToCompute*nParametersPerInstance_ << R"(;  // = nInstancesToCompute*)" << nParametersPerInstance_<< R"(;
const int nElementLengths = )" << (nInstancesToComputePerFiber_-1)*nFibersToCompute << R"(;  // = nElementsOnFiber*nFibersToCompute;
const int nFiringEvents = )" << gpuFiringEventsNRows_*gpuFiringEventsNColumns_ << R"(;  // = firingEventsNRows*firingEventsNColumns;
const int nFrequencyJitter = )" << nFibersToCompute*gpuFrequencyJitterNColumns_ << R"(;  // = nFibersToCompute*frequencyJitterNColumns;
const int nAlgebraicsForTransferIndices = )" << algebraicsForTransferIndices_.size() << R"(;
const int nAlgebraicsForTransfer = )" << nInstancesToCompute*algebraicsForTransferIndices_.size() << R"(;  // = nInstancesToCompute*nAlgebraicsForTransferIndices;
const int nStatesForTransferIndices = )" << statesForTransferIndices_.size() << R"(;
const int nStatesForTransfer = )" << nInstancesToCompute*statesForTransferIndices_.size() << R"(;  // = nInstancesToCompute*nStatesForTransferIndices;
)";

  if (optimizationType_ == "gpu")
  {
    sourceCode
      << "// the following code is generated by FastMonodomainSolver for offloading to GPU for "
      << 1 << " fibers, " << nInstancesToComputePerFiber_ << " instances per fiber, " << nInstancesToCompute_ << " instances in total.\n"
      << R"(

// global variables to be stored on the target device
#pragma omp declare target
double states[nStatesTotal];             // including state 0 which is stored in vmValues)";
  if (!algebraicsForTransferIndices_.empty())
  {
    sourceCode << R"(
int algebraicsForTransferIndices[nAlgebraicsForTransferIndices];
)";
  }
  sourceCode << R"(
int statesForTransferIndices[nStatesForTransferIndices];
char firingEvents[nFiringEvents];
double setSpecificStatesFrequencyJitter[nFrequencyJitter];
char fiberIsCurrentlyStimulated[nFibersToCompute];
int motorUnitNo[nFibersToCompute];
int fiberStimulationPointIndex[nFibersToCompute];
double lastStimulationCheckTime[nFibersToCompute];
double setSpecificStatesCallFrequency[nFibersToCompute];
double setSpecificStatesRepeatAfterFirstCall[nFibersToCompute];
double setSpecificStatesCallEnableBegin[nFibersToCompute];
double currentJitter[nFibersToCompute];
int jitterIndex[nFibersToCompute];

double vmValues[nInstancesToCompute];
double rates[nStatesTotal];
double intermediateRates[nStatesTotal];
double algebraics[nAlgebraicsTotal];
double intermediateAlgebraics[nAlgebraicsTotal];
#pragma omp end declare target

#ifdef __cplusplus
extern "C"
#endif
void initializeArrays(const double *statesOneInstance, const int *algebraicsForTransferIndicesParameter, const int *statesForTransferIndicesParameter,
                      const char *firingEventsParameter, const double *setSpecificStatesFrequencyJitterParameter, const int *motorUnitNoParameter,
                      const int *fiberStimulationPointIndexParameter, const double *lastStimulationCheckTimeParameter,
                      const double *setSpecificStatesCallFrequencyParameter, const double *setSpecificStatesRepeatAfterFirstCallParameter,
                      const double *setSpecificStatesCallEnableBeginParameter)
{
  // copy given values to variables on target
  for (int fiberNo = 0; fiberNo < nFibersToCompute; fiberNo++)
  {
    for (int instanceNo = 0; instanceNo < nInstancesPerFiber; instanceNo++)
    {
      int instanceToComputeNo = fiberNo*nInstancesPerFiber + instanceNo;

      // The entries in states[0] to states[1*nInstancesToCompute - 1] are not used.
      // State zero is stored in vmValues instead.
      for (int stateNo = 0; stateNo < nStates; stateNo++)
      {
        states[stateNo*nInstancesToCompute + instanceToComputeNo] = statesOneInstance[stateNo];
        rates[stateNo*nInstancesToCompute + instanceToComputeNo] = 0;
        intermediateRates[stateNo*nInstancesToCompute + instanceToComputeNo] = 0;
      }
    }
  }
)";
  if (!algebraicsForTransferIndices_.empty())
  {
    sourceCode << R"(
  for (int i = 0; i < nAlgebraicsForTransferIndices; i++)
    algebraicsForTransferIndices[i] = algebraicsForTransferIndicesParameter[i];
)";
  }
  sourceCode << R"(

  for (int i = 0; i < nStatesForTransferIndices; i++)
    statesForTransferIndices[i] = statesForTransferIndicesParameter[i];

  for (int i = 0; i < nFiringEvents; i++)
    firingEvents[i] = firingEventsParameter[i];

  for (int i = 0; i < nFrequencyJitter; i++)
    setSpecificStatesFrequencyJitter[i] = setSpecificStatesFrequencyJitterParameter[i];

  for (int fiberNo = 0; fiberNo < nFibersToCompute; fiberNo++)
  {
    motorUnitNo[fiberNo] = motorUnitNoParameter[fiberNo];
    fiberStimulationPointIndex[fiberNo] = fiberStimulationPointIndexParameter[fiberNo];
    lastStimulationCheckTime[fiberNo] = lastStimulationCheckTimeParameter[fiberNo];
    setSpecificStatesCallFrequency[fiberNo] = setSpecificStatesCallFrequencyParameter[fiberNo];
    setSpecificStatesRepeatAfterFirstCall[fiberNo] = setSpecificStatesRepeatAfterFirstCallParameter[fiberNo];
    setSpecificStatesCallEnableBegin[fiberNo] = setSpecificStatesCallEnableBeginParameter[fiberNo];
  }

  // set variables to zero
  for (int fiberNo = 0; fiberNo < nFibersToCompute; fiberNo++)
  {
    fiberIsCurrentlyStimulated[fiberNo] = 0;
    currentJitter[fiberNo] = 0;
    jitterIndex[fiberNo] = 0;
  }

  // initialize vmValues
  const double state0 = statesOneInstance[0];
  for (int instanceToComputeNo = 0; instanceToComputeNo < nInstancesToCompute; instanceToComputeNo++)
  {
    vmValues[instanceToComputeNo] = state0;
  }

)";
  if (!algebraicsForTransferIndices_.empty())
  {
    sourceCode << R"(
  // map values to target
  #pragma omp target update to(states[:nStatesTotal], algebraicsForTransferIndices[:nAlgebraicsForTransferIndices], statesForTransferIndices[:nStatesForTransferIndices], \
    firingEvents[:nFiringEvents], setSpecificStatesFrequencyJitter[:nFrequencyJitter], \
    motorUnitNo[:nFibersToCompute], fiberStimulationPointIndex[:nFibersToCompute], \
    lastStimulationCheckTime[:nFibersToCompute], setSpecificStatesCallFrequency[:nFibersToCompute], \
    setSpecificStatesRepeatAfterFirstCall[:nFibersToCompute], setSpecificStatesCallEnableBegin[:nFibersToCompute], \
    currentJitter[:nFibersToCompute], jitterIndex[:nFibersToCompute], vmValues[:nInstancesToCompute]))";
  }
  else
  {
    sourceCode << R"(
  // map values to target
  #pragma omp target update to(states[:nStatesTotal], statesForTransferIndices[:nStatesForTransferIndices], \
    firingEvents[:nFiringEvents], setSpecificStatesFrequencyJitter[:nFrequencyJitter], \
    motorUnitNo[:nFibersToCompute], fiberStimulationPointIndex[:nFibersToCompute], \
    lastStimulationCheckTime[:nFibersToCompute], setSpecificStatesCallFrequency[:nFibersToCompute], \
    setSpecificStatesRepeatAfterFirstCall[:nFibersToCompute], setSpecificStatesCallEnableBegin[:nFibersToCompute], \
    currentJitter[:nFibersToCompute], jitterIndex[:nFibersToCompute], vmValues[:nInstancesToCompute]))";
  }
  sourceCode << R"(
}

)";

  }

  sourceCode
    << "// compute the total monodomain equation\n"
    << "#ifdef __cplusplus\n" << "extern \"C\"\n" << "#endif"
    << R"(
void computeMonodomain(const double *parameters,
                       double *algebraicsForTransfer, double *statesForTransfer, const double *elementLengths,
                       double startTime, double timeStepWidthSplitting, int nTimeStepsSplitting, double dt0D, int nTimeSteps0D, double dt1D, int nTimeSteps1D,
                       double prefactor, double valueForStimulatedPoint)
{
)";

  if (optimizationType_ == "gpu")
    sourceCode << R"(

  // map data to and from GPU
  #pragma omp target data map(to: parameters[:nParametersTotal], elementLengths[:nElementLengths]) \
       map(from: )";
  if (!algebraicsForTransferIndices_.empty())
  {
    sourceCode << R"(algebraicsForTransfer[:nAlgebraicsForTransfer], )";
  }
  sourceCode << R"(statesForTransfer[:nStatesForTransfer])
  {
)";
   sourceCode << R"(
  //const int nAlgebraics = )" << nAlgebraics << R"(;

  // loop over splitting time steps
  #pragma omp target teams
  for (int timeStepNo = 0; timeStepNo < nTimeStepsSplitting; timeStepNo++)
  {
    // perform Strang splitting
    double currentTimeSplitting = startTime + timeStepNo * timeStepWidthSplitting;

    // compute midTimeSplitting once per step to reuse it. [currentTime, midTimeSplitting=currentTime+0.5*timeStepWidth, currentTime+timeStepWidth]
    double midTimeSplitting = currentTimeSplitting + 0.5 * timeStepWidthSplitting;
    bool storeAlgebraicsForTransferSplitting = false;   // do not store the computed algebraics values in the algebraicsForTransfer vector for communication, because this is the first 0D computation and they will be changed in the second 0D computation

    // perform Strang splitting:
    // 0D: [currentTimeSplitting, currentTimeSplitting + dt0D*nTimeSteps0D]
    // 1D: [currentTimeSplitting, currentTimeSplitting + dt1D*nTimeSteps1D]
    // 0D: [midTimeSplitting,     midTimeSplitting + dt0D*nTimeSteps0D]

    // advance 0D in [currentTimeSplitting, currentTimeSplitting + dt0D*nTimeSteps0D]
    // ------------------------------------------------------------

    // loop over fibers that will be computed on this rank)";
  if (optimizationType_ == "openmp")
    sourceCode << "\n    #pragma omp parallel for";
  else if (optimizationType_ == "gpu")
    sourceCode << "\n    #pragma omp distribute parallel for simd collapse(2)";  // teams distribute
  sourceCode << R"(
    for (int fiberNo = 0; fiberNo < nFibersToCompute; fiberNo++)
    {
      // loop over instances to compute here
      for (int instanceNo = 0; instanceNo < nInstancesPerFiber; instanceNo++)
      {
        int instanceToComputeNo = fiberNo*nInstancesPerFiber + instanceNo;    // index of instance over all fibers

        // determine if current point is at center of fiber
        int fiberCenterIndex = fiberStimulationPointIndex[fiberNo];
        bool currentPointIsInCenter = fabs(fiberCenterIndex - instanceNo) < 4;

        // loop over 0D timesteps
        for (int timeStepNo = 0; timeStepNo < nTimeSteps0D; timeStepNo++)
        {
          double currentTime = currentTimeSplitting + timeStepNo * dt0D;

          // determine if fiber gets stimulated
          // check if current point will be stimulated
          bool stimulateCurrentPoint = false;
          if (currentPointIsInCenter)
          {
            // check if time has come to call setSpecificStates
            bool checkStimulation = false;

            if (currentTime >= lastStimulationCheckTime[fiberNo] + 1./(setSpecificStatesCallFrequency[fiberNo]+currentJitter[fiberNo])
                && currentTime >= setSpecificStatesCallEnableBegin[fiberNo]-1e-13)
            {
              checkStimulation = true;

              // if current stimulation is over
              if (setSpecificStatesRepeatAfterFirstCall[fiberNo] != 0
                  && currentTime - (lastStimulationCheckTime[fiberNo] + 1./(setSpecificStatesCallFrequency[fiberNo] + currentJitter[fiberNo])) > setSpecificStatesRepeatAfterFirstCall[fiberNo])
              {
                // advance time of last call to specificStates
                lastStimulationCheckTime[fiberNo] += 1./(setSpecificStatesCallFrequency[fiberNo] + currentJitter[fiberNo]);

                // compute new jitter value
                double jitterFactor = 0.0;
                if (frequencyJitterNColumns > 0)
                  jitterFactor = setSpecificStatesFrequencyJitter[fiberNo*frequencyJitterNColumns + jitterIndex[fiberNo] % frequencyJitterNColumns];
                currentJitter[fiberNo] = jitterFactor * setSpecificStatesCallFrequency[fiberNo];

                jitterIndex[fiberNo]++;

                checkStimulation = false;
              }
            }

            // instead of calling setSpecificStates, directly determine whether to stimulate from the firingEvents file
            int firingEventsTimeStepNo = int(currentTime * setSpecificStatesCallFrequency[fiberNo] + 0.5);
            int firingEventsIndex = (firingEventsTimeStepNo % firingEventsNRows)*firingEventsNColumns + (motorUnitNo[fiberNo] % firingEventsNColumns);
            // firingEvents_[timeStepNo*nMotorUnits + motorUnitNo[fiberNo]]

            stimulateCurrentPoint = checkStimulation && firingEvents[firingEventsIndex];
            fiberIsCurrentlyStimulated[fiberNo] = stimulateCurrentPoint? 1: 0;

            // output to console
            if (stimulateCurrentPoint && fiberCenterIndex == instanceNo)
            {
              if (omp_is_initial_device())
                printf("t: %f, stimulate fiber %d (local no.), MU %d (computation on CPU)\n", currentTime, fiberNo, motorUnitNo[fiberNo]);
              else
                printf("t: %f, stimulate fiber %d (local no.), MU %d (computation on GPU)\n", currentTime, fiberNo, motorUnitNo[fiberNo]);
            }
          }
          const bool storeAlgebraicsForTransfer = false;

          )" << mainCode << R"(        }  // loop over 0D timesteps
      }  // loop over instances)";

  if (optimizationType_ == "gpu")
    sourceCode << R"(
    }  // loop over fibers
)";

    // depending on DiffusionTimeSteppingScheme either do Implicit Euler or Crank-Nicolson
    bool useImplicitEuler = std::is_same<DiffusionTimeSteppingScheme,
                              TimeSteppingScheme::ImplicitEuler<typename DiffusionTimeSteppingScheme::DiscretizableInTime>
                            >::value;

    sourceCode << R"(
    // advance 1D in [currentTimeSplitting, currentTimeSplitting + dt1D*nTimeSteps1D]
    // ------------------------------------------------------------

    // Implicit Euler step:
    // (K - 1/dt*M) u^{n+1} = -1/dt*M u^{n})
    // Crank-Nicolson step:
    // (1/2*K - 1/dt*M) u^{n+1} = (-1/2*K -1/dt*M) u^{n})

    // stencil K: 1/h*[_-1_  1  ]*prefactor
    // stencil M:   h*[_1/3_ 1/6]
    const double dt = dt1D;
)";
  if (optimizationType_ == "gpu")
    sourceCode << R"(
    // loop over fibers
    #pragma omp distribute parallel for
    for (int fiberNo = 0; fiberNo < nFibersToCompute; fiberNo++)
    {
)";
  sourceCode << R"(
      const int nValues = nInstancesPerFiber;

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

      // helper buffers c', d' for Thomas algorithm
      double cIntermediate[nInstancesPerFiber-1];
      double dIntermediate[nInstancesPerFiber];

      // perform forward substitution
      // loop over entries / rows of matrices, this is equal to the instances of the current fiber
      for (int valueNo = 0; valueNo < nValues; valueNo++)
      {
        int instanceToComputeNo = fiberNo*nInstancesPerFiber + valueNo;

        double a = 0;
        double b = 0;
        double c = 0;
        double d = 0;

        double u_previous = 0;
        double u_center = 0;
        double u_next = 0;

        u_center = vmValues[instanceToComputeNo];  // state 0 of the current instance

        // contribution from left element
        if (valueNo > 0)
        {
          u_previous = vmValues[instanceToComputeNo - 1];  // state 0 of the left instance

          // stencil K: 1/h*[1   _-1_ ]*prefactor
          // stencil M:   h*[1/6 _1/3_]

          double h_left = elementLengths[fiberNo*nElementsOnFiber + valueNo-1];
          double k_left = 1./h_left*(1) * prefactor;
          double m_left = h_left*1./6;
)";
  if (useImplicitEuler)
  {
    sourceCode << R"(
          a = (k_left - 1/dt*m_left);   // formula for implicit Euler
)";
  }
  else  // Crank-Nicolson
  {
    sourceCode << R"(
          a = (k_left/2. - 1/dt*m_left);   // formula for Crank-Nicolson
)";
  }

  sourceCode << R"(
          double k_right = 1./h_left*(-1) * prefactor;
          double m_right = h_left*1./3;
)";
  if (useImplicitEuler)
  {
    sourceCode << R"(
          b += (k_right - 1/dt*m_right);     // formula for implicit Euler
          d += (-1/dt*m_left) * u_previous + (-1/dt*m_right) * u_center;
)";
  }
  else  // Crank-Nicolson
  {
    sourceCode << R"(
          b += (k_right/2. - 1/dt*m_right);   // formula for Crank-Nicolson
          d += (-k_left/2. - 1/dt*m_left) * u_previous + (-k_right/2. - 1/dt*m_right) * u_center;
)";
  }
  sourceCode << R"(
        }

        // contribution from right element
        if (valueNo < nValues-1)
        {
          u_next = vmValues[instanceToComputeNo + 1];  // state 0 of the right instance

          // stencil K: 1/h*[_-1_  1  ]*prefactor
          // stencil M:   h*[_1/3_ 1/6]

          double h_right = elementLengths[fiberNo*nElementsOnFiber + valueNo];
          double k_right = 1./h_right*(1) * prefactor;
          double m_right = h_right*1./6;
)";
  if (useImplicitEuler)
  {
    sourceCode << R"(
          c = (k_right - 1/dt*m_right);     // formula for implicit Euler
)";
  }
  else  // Crank-Nicolson
  {
    sourceCode << R"(
          c = (k_right/2. - 1/dt*m_right);   // formula for Crank-Nicolson
)";
  }
  sourceCode << R"(

          double k_left = 1./h_right*(-1) * prefactor;
          double m_left = h_right*1./3;
)";

  if (useImplicitEuler)
  {
    sourceCode << R"(
          b += (k_left - 1/dt*m_left);
          d += (-1/dt*m_left) * u_center + (-1/dt*m_right) * u_next;     // formula for implicit Euler
)";
  }
  else  // Crank-Nicolson
  {
    sourceCode << R"(
          b += (k_left/2. - 1/dt*m_left);
          d += (-k_left/2. - 1/dt*m_left) * u_center + (-k_right/2. - 1/dt*m_right) * u_next;   // formula for Crank-Nicolson
)";
  }
  sourceCode << R"(
        }

        if (valueNo == 0)
        {
          // c'_0 = c_0 / b_0
          cIntermediate[valueNo] = c / b;

          // d'_0 = d_0 / b_0
          dIntermediate[valueNo] = d / b;
        }
        else
        {
          if (valueNo != nValues-1)
          {
            // c'_i = c_i / (b_i - c'_{i-1}*a_i)
            cIntermediate[valueNo] = c / (b - cIntermediate[valueNo-1]*a);
          }

          // d'_i = (d_i - d'_{i-1}*a_i) / (b_i - c'_{i-1}*a_i)
          dIntermediate[valueNo] = (d - dIntermediate[valueNo-1]*a) / (b - cIntermediate[valueNo-1]*a);
        }
      }

      // perform backward substitution
      // x_n = d'_n
      vmValues[nValues-1] = dIntermediate[nValues-1];  // state 0 of the point (nValues-1)

      double previousValue = dIntermediate[nValues-1];

      // loop over entries / rows of matrices
      for (int valueNo = nValues-2; valueNo >= 0; valueNo--)
      {
        int instanceToComputeNo = fiberNo*nInstancesPerFiber + valueNo;

        // x_i = d'_i - c'_i * x_{i+1}
        double resultValue = dIntermediate[valueNo] - cIntermediate[valueNo] * previousValue;
        vmValues[instanceToComputeNo] = resultValue;

        previousValue = resultValue;
      }
)";
  if (optimizationType_ == "gpu")
    sourceCode << "    }";
  sourceCode << R"(

    // advance 0D in [midTimeSplitting,     midTimeSplitting + dt0D*nTimeSteps0D]
    // ------------------------------------------------------------
    // in the last timestep, store the computed algebraics values in the algebraicsForTransfer vector for communication
    storeAlgebraicsForTransferSplitting = timeStepNo == nTimeStepsSplitting-1;
)";
  if (optimizationType_ == "gpu")
    sourceCode << R"(
    // loop over fibers that will be computed on this rank

    #pragma omp distribute parallel for simd collapse(2)
    for (int fiberNo = 0; fiberNo < nFibersToCompute; fiberNo++)
    {)";
  sourceCode << R"(
      // loop over instances to compute here
      for (int instanceNo = 0; instanceNo < nInstancesPerFiber; instanceNo++)
      {
        int instanceToComputeNo = fiberNo*nInstancesPerFiber + instanceNo;    // index of instance over all fibers

        // determine if current point is at center of fiber
        int fiberCenterIndex = fiberStimulationPointIndex[fiberNo];
        bool currentPointIsInCenter = fabs(fiberCenterIndex - instanceNo) < 4;

        // loop over 0D timesteps
        for (int timeStepNo = 0; timeStepNo < nTimeSteps0D; timeStepNo++)
        {
          double currentTime = midTimeSplitting + timeStepNo * dt0D;

          // determine if fiber gets stimulated
          // check if current point will be stimulated
          bool stimulateCurrentPoint = false;
          if (currentPointIsInCenter)
          {
            // check if time has come to call setSpecificStates
            bool checkStimulation = false;

            if (currentTime >= lastStimulationCheckTime[fiberNo] + 1./(setSpecificStatesCallFrequency[fiberNo]+currentJitter[fiberNo])
                && currentTime >= setSpecificStatesCallEnableBegin[fiberNo]-1e-13)
            {
              checkStimulation = true;

              // if current stimulation is over
              if (setSpecificStatesRepeatAfterFirstCall[fiberNo] != 0
                  && currentTime - (lastStimulationCheckTime[fiberNo] + 1./(setSpecificStatesCallFrequency[fiberNo] + currentJitter[fiberNo])) > setSpecificStatesRepeatAfterFirstCall[fiberNo])
              {
                // advance time of last call to specificStates
                lastStimulationCheckTime[fiberNo] += 1./(setSpecificStatesCallFrequency[fiberNo] + currentJitter[fiberNo]);

                // compute new jitter value
                double jitterFactor = 0.0;
                if (frequencyJitterNColumns > 0)
                  jitterFactor = setSpecificStatesFrequencyJitter[fiberNo*frequencyJitterNColumns + jitterIndex[fiberNo] % frequencyJitterNColumns];
                currentJitter[fiberNo] = jitterFactor * setSpecificStatesCallFrequency[fiberNo];

                jitterIndex[fiberNo]++;

                checkStimulation = false;
              }
            }

            // instead of calling setSpecificStates, directly determine whether to stimulate from the firingEvents file
            int firingEventsTimeStepNo = int(currentTime * setSpecificStatesCallFrequency[fiberNo] + 0.5);
            int firingEventsIndex = (firingEventsTimeStepNo % firingEventsNRows)*firingEventsNColumns + (motorUnitNo[fiberNo] % firingEventsNColumns);
            // firingEvents_[timeStepNo*nMotorUnits + motorUnitNo[fiberNo]]

            stimulateCurrentPoint = checkStimulation && firingEvents[firingEventsIndex];
            fiberIsCurrentlyStimulated[fiberNo] = stimulateCurrentPoint? 1: 0;
          }
          const bool storeAlgebraicsForTransfer = storeAlgebraicsForTransferSplitting && timeStepNo == nTimeSteps0D-1;

          )" << mainCode << R"(

        }  // loop over 0D timesteps
      }  // loop over instances
    }  // loop over fibers
  } // loop over splitting timesteps
)";
  if (optimizationType_ == "gpu")
    sourceCode << R"(
  } // end pragma omp target
)";
  sourceCode << R"(
  // map back from GPU to host
  //#pragma omp target update from()";
  if (!algebraicsForTransferIndices_.empty())
  {
    sourceCode << R"(algebraicsForTransfer[:nAlgebraicsForTransfer], )";
  }
  sourceCode << R"(statesForTransfer[:nStatesForTransfer])

}
)";
    
  // write to source file
  std::ofstream sourceCodeFile;
  OutputWriter::Generic::openFile(sourceCodeFile, outputFilename);
  if (!sourceCodeFile.is_open())
  {
    LOG(FATAL) << "Could not write to file \"" << outputFilename << "\".";
  }
  else
  {
    std::string fileContents = sourceCode.str();
    sourceCodeFile << fileContents;
    sourceCodeFile.close();
  }
}

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
initializeValuesOnGpu()
{
  // check that arrays have the correct sizes
  /*
  // size constants
  const int nElementsOnFiber = nInstancesToComputePerFiber_-1;
  const int nFibersToCompute = 1;
  const int nParametersTotal = nInstancesToCompute_*nParametersPerInstance_;
  const int nElementLengths = nElementsOnFiber*nFibersToCompute;
  const int nFiringEvents = gpuFiringEventsNRows_*gpuFiringEventsNColumns_;
  const int nFrequencyJitter = nFibersToCompute*gpuFrequencyJitterNColumns_;
  const int nAlgebraicsForTransferIndices = algebraicsForTransferIndices_.size();
  const int nAlgebraicsForTransfer = nInstancesToCompute_*nAlgebraicsForTransferIndices;
  const int nStatesForTransferIndices = statesForTransferIndices_.size();
  const int nStatesForTransfer = nInstancesToCompute_*nStatesForTransferIndices;

  if (gpuVmValues_.size() != nInstancesToCompute_)
    LOG(FATAL) << "gpuVmValues_.size() = " << gpuVmValues_.size() << " does not match assumed size " << nInstancesToCompute_;
  if (gpuFiberIsCurrentlyStimulated_.size() != nFibersToCompute)
    LOG(FATAL) << "gpuFiberIsCurrentlyStimulated_.size() = " << gpuFiberIsCurrentlyStimulated_.size() << " does not match assumed size " << nFibersToCompute;
  if (gpuLastStimulationCheckTime_.size() != nFibersToCompute)
    LOG(FATAL) << "gpuLastStimulationCheckTime_.size() = " << gpuLastStimulationCheckTime_.size() << " does not match assumed size " << nFibersToCompute;
  if (gpuCurrentJitter_.size() != nFibersToCompute)
    LOG(FATAL) << "gpuCurrentJitter_.size() = " << gpuCurrentJitter_.size() << " does not match assumed size " << nFibersToCompute;
  if (gpuJitterIndex_.size() != nFibersToCompute)
    LOG(FATAL) << "gpuJitterIndex_.size() = " << gpuJitterIndex_.size() << " does not match assumed size " << nFibersToCompute;

  if (gpuParameters_.size() != nParametersTotal)
    LOG(FATAL) << "gpuParameters_.size() = " << gpuParameters_.size() << " does not match assumed size " << nParametersTotal << "x" << nInstancesToCompute_ << "x" << nParametersPerInstance_;
  if (algebraicsForTransferIndices_.size() != nAlgebraicsForTransferIndices)
    LOG(FATAL) << "algebraicsForTransferIndices_.size() = " << algebraicsForTransferIndices_.size() << " does not match assumed size " << nAlgebraicsForTransferIndices;
  if (statesForTransferIndices_.size() != nStatesForTransferIndices)
    LOG(FATAL) << "statesForTransferIndices_.size() = " << statesForTransferIndices_.size() << " does not match assumed size " << nStatesForTransferIndices;
  if (gpuElementLengths_.size() != nElementLengths)
    LOG(FATAL) << "gpuElementLengths_.size() = " << gpuElementLengths_.size() << " does not match assumed size " << nElementLengths << "=" << nElementsOnFiber << "x" << nFibersToCompute;
  if (gpuFiringEvents_.size() != nFiringEvents)
    LOG(FATAL) << "gpuFiringEvents_.size() = " << gpuFiringEvents_.size() << " does not match assumed size " << nFiringEvents << "=" << gpuFiringEventsNRows_ << "x" << gpuFiringEventsNColumns_;
  if (gpuSetSpecificStatesFrequencyJitter_.size() != nFrequencyJitter)
    LOG(FATAL) << "gpuSetSpecificStatesFrequencyJitter_.size() = " << gpuSetSpecificStatesFrequencyJitter_.size() << " does not match assumed size " << nFrequencyJitter << "=" << 1 << "x" << gpuFrequencyJitterNColumns_;
  if (gpuMotorUnitNo_.size() != nFibersToCompute)
    LOG(FATAL) << "gpuMotorUnitNo_.size() = " << gpuMotorUnitNo_.size() << " does not match assumed size " << nFibersToCompute;

  if (gpuMotorUnitNo_.size() != nFibersToCompute)
    LOG(FATAL) << "gpuMotorUnitNo_.size() = " << gpuMotorUnitNo_.size() << " does not match assumed size " << nFibersToCompute;
  if (gpuFiberStimulationPointIndex_.size() != nFibersToCompute)
    LOG(FATAL) << "gpuFiberStimulationPointIndex_.size() = " << gpuFiberStimulationPointIndex_.size() << " does not match assumed size " << nFibersToCompute;
  if (gpuSetSpecificStatesCallFrequency_.size() != nFibersToCompute)
    LOG(FATAL) << "gpuSetSpecificStatesCallFrequency_.size() = " << gpuSetSpecificStatesCallFrequency_.size() << " does not match assumed size " << nFibersToCompute;
  if (gpuSetSpecificStatesRepeatAfterFirstCall_.size() != nFibersToCompute)
    LOG(FATAL) << "gpuSetSpecificStatesRepeatAfterFirstCall_.size() = " << gpuSetSpecificStatesRepeatAfterFirstCall_.size() << " does not match assumed size " << nFibersToCompute;
  if (gpuSetSpecificStatesCallEnableBegin_.size() != nFibersToCompute)
    LOG(FATAL) << "gpuSetSpecificStatesCallEnableBegin_.size() = " << gpuSetSpecificStatesCallEnableBegin_.size() << " does not match assumed size " << nFibersToCompute;

  if (gpuAlgebraicsForTransfer_.size() != nAlgebraicsForTransfer)
    LOG(FATAL) << "gpuAlgebraicsForTransfer_.size() = " << gpuAlgebraicsForTransfer_.size() << " does not match assumed size " << nAlgebraicsForTransfer << "=" << nInstancesToCompute_ << "x" << nAlgebraicsForTransferIndices;
  if (gpuStatesForTransfer_.size() != nStatesForTransfer)
    LOG(FATAL) << "gpuStatesForTransfer_.size() = " << gpuStatesForTransfer_.size() << " does not match assumed size " << nStatesForTransfer << "=" << nInstancesToCompute_ << "x" << nStatesForTransferIndices;
*/
  CellmlAdapterType &cellmlAdapter = nestedSolvers_.instancesLocal()[0].timeStepping1().instancesLocal()[0].discretizableInTime();
  CellmlSourceCodeGenerator &cellmlSourceCodeGenerator = cellmlAdapter.cellmlSourceCodeGenerator();
  const std::vector<double> &statesInitialValues = cellmlSourceCodeGenerator.statesInitialValues();

  // upload all values to GPU
  initializeArrays_(statesInitialValues.data(), algebraicsForTransferIndices_.data(), statesForTransferIndices_.data(),
                   gpuFiringEvents_.data(), gpuSetSpecificStatesFrequencyJitter_.data(), gpuMotorUnitNo_.data(),
                   gpuFiberStimulationPointIndex_.data(), gpuLastStimulationCheckTime_.data(),
                   gpuSetSpecificStatesCallFrequency_.data(), gpuSetSpecificStatesRepeatAfterFirstCall_.data(),
                   gpuSetSpecificStatesCallEnableBegin_.data());
}

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
computeMonodomainGpu()
{
  LOG(TRACE) << "computeMonodomainGpu";

  // initialize data vector

  // fetch timestep widths and total time span
  std::vector<typename NestedSolversType::TimeSteppingSchemeType> &instances = nestedSolvers_.instancesLocal();

  TimeSteppingScheme::Heun<CellmlAdapterType> &heun = instances[0].timeStepping1().instancesLocal()[0];

  DiffusionTimeSteppingScheme &implicitEuler = instances[0].timeStepping2().instancesLocal()[0];
  double prefactor = implicitEuler.discretizableInTime().data().context().getPythonConfig().getOptionDouble("prefactor", 1.0);

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

  if (fiberData_.empty())
  {
    LOG(DEBUG) << "In computeMonodomain(" << startTime << "," << timeStepWidthSplitting
      << ") fiberData_ is empty";
    LOG(DEBUG) << "This means there is no fiber to compute on this rank, they were all send to another rank for the computation. Skip computation.";
    return;
  }

  if (!computeMonodomain_)
    LOG(FATAL) << "Could not call computeMonodomain function.";

  if (nInstancesToCompute_ == 0)
    return;

  // copy all states values from fiberData to the gpuVmValues vector which is then used in the generated code
  for (int fiberDataNo = 0; fiberDataNo < 1; fiberDataNo++)
  {
    for (int instanceNo = 0; instanceNo < nInstancesToComputePerFiber_; instanceNo++)
    {
      const int instanceNoTotal = fiberDataNo*nInstancesToComputePerFiber_ + instanceNo;

      // set state 0
      gpuVmValues_[instanceNoTotal] = fiberData_[fiberDataNo].vmValues[instanceNo];
    }
  }

  // call the compiled function
  computeMonodomain_(gpuParameters_.data(),
                     gpuAlgebraicsForTransfer_.data(), gpuStatesForTransfer_.data(), gpuElementLengths_.data(),
                     startTime, timeStepWidthSplitting, nTimeStepsSplitting_, dt0D, nTimeSteps0D, dt1D, nTimeSteps1D,
                     prefactor, valueForStimulatedPoint_);

  // copy the resulting values back to fiberData_
  for (int fiberDataNo = 0; fiberDataNo < 1; fiberDataNo++)
  {
    for (int instanceNo = 0; instanceNo < nInstancesToComputePerFiber_; instanceNo++)
    {
      const int instanceNoTotal = fiberDataNo*nInstancesToComputePerFiber_ + instanceNo;

      if (statesForTransferIndices_.size() == 0)
        LOG(FATAL) << "The implementation of the FastMonodomainSolver needs the transmembrane potential to be transferred.";

      // the first state to transfer is the Vm value
      fiberData_[fiberDataNo].vmValues[instanceNo] = gpuStatesForTransfer_[0*nInstancesToCompute_ + instanceNoTotal];

      // loop over further states to transfer
      int furtherDataIndex = 0;
      for (int i = 1; i < statesForTransferIndices_.size(); i++, furtherDataIndex++)
      {
        fiberData_[fiberDataNo].furtherStatesAndAlgebraicsValues[furtherDataIndex*nInstancesToComputePerFiber_ + instanceNo]
          = gpuStatesForTransfer_[i*nInstancesToCompute_ + instanceNoTotal];
      }

      // loop over algebraics to transfer
      for (int i = 0; i < algebraicsForTransferIndices_.size(); i++, furtherDataIndex++)
      {
        fiberData_[fiberDataNo].furtherStatesAndAlgebraicsValues[furtherDataIndex*nInstancesToComputePerFiber_ + instanceNo]
          = gpuAlgebraicsForTransfer_[i*nInstancesToCompute_ + instanceNoTotal];
      }
    }
    LOG(DEBUG) << "states and algebraics for transfer at fiberDataNo=" << fiberDataNo << ": " << fiberData_[fiberDataNo].furtherStatesAndAlgebraicsValues;
    LOG(DEBUG) << "size: " << fiberData_[fiberDataNo].furtherStatesAndAlgebraicsValues.size() << ", nInstancesToComputePerFiber_: " << nInstancesToComputePerFiber_;
  }

  currentTime_ = instances[0].endTime();
}
