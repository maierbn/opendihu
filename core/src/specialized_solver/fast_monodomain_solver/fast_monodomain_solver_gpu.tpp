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
  cellmlSourceCodeGenerator.generateSourceFastMonodomainGpu(approximateExponentialFunction,
                                                            nFibersToCompute_, nInstancesToComputePerFiber_, nParametersPerInstance_,
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

  //LOG(ERROR) << "generateMonodomainSolverGpuSource commented out";
  generateMonodomainSolverGpuSource(sourceToCompileFilename, headerCode, mainCode);

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

  computeMonodomain_ = (void (*)(FiberData *fiberData, double *states, const double *parameters,
                             double *algebraicsForTransfer, const int *algebraicsForTransferIndices, int nAlgebraicsForTransferIndices,
                             const double *elementLengths, char *firingEvents, int firingEventsNRows, int firingEventsNColumns,
                             double *setSpecificStatesFrequencyJitter, int frequencyJitterNColumns, char *fiberIsCurrentlyStimulated,
                             double startTime, double timeStepWidthSplitting, int nTimeStepsSplitting, double dt0D, int nTimeSteps0D, double dt1D, int nTimeSteps1D,
                             double prefactor, double valueForStimulatedPoint)) dlsym(handle, "computeMonodomain");

  LOG(DEBUG) << "computeMonodomain: " << (computeMonodomain_==nullptr? "no" : "yes");

  if (computeMonodomain_ == nullptr)
  {
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
  
  sourceCode << "// the following code is generated by FastMonodomainSolver for offloading to GPU for "
    << nFibersToCompute_ << " fibers, " << nInstancesToComputePerFiber_ << " instances per fiber, " << nInstancesToCompute_ << " instances in total.";

  // define the struct
  sourceCode << R"(
typedef unsigned long long global_no_t;

// struct with which data is transferred from the calling code
struct FiberData
{
  std::vector<double> elementLengths;   //< lengths of the 1D elements
  std::vector<double> vmValues;         //< values of Vm
  std::vector<double> furtherStatesAndIntermediatesValues;    //< all data to be transferred back to the fibers, apart from vmValues, corresponding to statesForTransferIndices_ and algebraicsForTransferIndices_ (array of struct memory layout)
  int valuesLength;                     //< number of vmValues
  global_no_t valuesOffset;             //< number of vmValues in previous entries in fiberData_

  int fiberNoGlobal;                    //< fiberNo as given in settings (value of additionalArgument)
  int motorUnitNo;                      //< motor unit no.
  int fiberStimulationPointIndex;       //< index of the point on the fiber where to stimulate, i.e. position of the neuromuscular junction, if at center, it is equal to (int)(fiberData_[fiberDataNo].valuesLength / 2)

  double lastStimulationCheckTime;      //< last time the fiber was checked for stimulation
  double setSpecificStatesCallFrequency;        //< value of option with the same name in the python settings
  std::vector<double> setSpecificStatesFrequencyJitter;      //< value of option with the same name in the python settings
  double setSpecificStatesRepeatAfterFirstCall; //< how long in ms the prescribed value should be set
  double setSpecificStatesCallEnableBegin;      //< value of option with the same name in the python settings

  double currentJitter;                         //< current absolute value of jitter to add to setSpecificStatesCallFrequency
  int jitterIndex;                              //< index of the vector in setSpecificStatesFrequencyJitter which is the current value to use
  bool currentlyStimulating;                    //< if a stimulation is in progress at the current time
};

)";

  // define compute0D which computes one Heun step
  sourceCode
    << "// determine if an instance is stimulated on a fiber at the current point in time\n"
    << "#ifdef __cplusplus\n" << "extern \"C\"\n" << "#endif" << std::endl
    << R"(bool isCurrentPointStimulated(FiberData *fiberData, char *firingEvents, int firingEventsNRows, int firingEventsNColumns,
                              double *setSpecificStatesFrequencyJitter, int frequencyJitterNColumns, char *fiberIsCurrentlyStimulated,
                              int fiberDataNo, double currentTime, bool currentPointIsInCenter)
{
  FiberData &fiberDataCurrentPoint = fiberData[fiberDataNo];

  // there is a parallel piece of code to this one in CellmlAdapter<>::checkCallbackStates, cellml/03_cellml_adapter.tpp

  // get time from with testing for stimulation is enabled
  const double lastStimulationCheckTime                 = fiberDataCurrentPoint.lastStimulationCheckTime;

  const double setSpecificStatesCallFrequency           = fiberDataCurrentPoint.setSpecificStatesCallFrequency;
  const double setSpecificStatesRepeatAfterFirstCall    = fiberDataCurrentPoint.setSpecificStatesRepeatAfterFirstCall;
  const double setSpecificStatesCallEnableBegin         = fiberDataCurrentPoint.setSpecificStatesCallEnableBegin;

  const int &motorUnitNo                                = fiberDataCurrentPoint.motorUnitNo;
  int &jitterIndex                                      = fiberDataCurrentPoint.jitterIndex;
  double &currentJitter                                 = fiberDataCurrentPoint.currentJitter;

  // check if time has come to call setSpecificStates
  bool checkStimulation = false;

  if (currentTime >= lastStimulationCheckTime + 1./(setSpecificStatesCallFrequency+currentJitter)
      && currentTime >= setSpecificStatesCallEnableBegin-1e-13)
  {
    checkStimulation = true;

    // if current stimulation is over
    if (setSpecificStatesRepeatAfterFirstCall != 0
        && currentTime - (lastStimulationCheckTime + 1./(setSpecificStatesCallFrequency+currentJitter)) > setSpecificStatesRepeatAfterFirstCall)
    {
      // advance time of last call to specificStates
      fiberDataCurrentPoint.lastStimulationCheckTime += 1./(setSpecificStatesCallFrequency+currentJitter);

      // compute new jitter value
      double jitterFactor = 0.0;
      if (frequencyJitterNColumns > 0)
        jitterFactor = setSpecificStatesFrequencyJitter[fiberDataNo*frequencyJitterNColumns + jitterIndex % frequencyJitterNColumns];
      currentJitter = jitterFactor * setSpecificStatesCallFrequency;

      jitterIndex++;

      checkStimulation = false;
    }
  }

  // instead of calling setSpecificStates, directly determine whether to stimulate from the firingEvents file
  int firingEventsTimeStepNo = round(currentTime * setSpecificStatesCallFrequency);
  int firingEventsIndex = (firingEventsTimeStepNo % firingEventsNRows)*firingEventsNColumns + (motorUnitNo % firingEventsNColumns);
  // firingEvents_[timeStepNo*nMotorUnits + motorUnitNo]

  bool stimulate = checkStimulation && firingEvents[firingEventsIndex];

  if (stimulate && currentPointIsInCenter)
  {
    // if this is the first point in time of the current stimulation, log stimulation time
    if (!fiberIsCurrentlyStimulated[fiberDataNo])
    {
      fiberIsCurrentlyStimulated[fiberDataNo] = 1;
    }
  }

  if (!stimulate)
  {
    fiberIsCurrentlyStimulated[fiberDataNo] = 0;
  }

  bool stimulateCurrentPoint = stimulate && currentPointIsInCenter;
  return stimulateCurrentPoint;
}

)";

  if (optimizationType_ == "gpu")
    sourceCode << "#pragma omp end declare target\n";

  sourceCode
    << "// compute the total monodomain equation\n"
    << "#ifdef __cplusplus\n" << "extern \"C\"\n" << "#endif"
    << R"(
void computeMonodomain(FiberData *fiberData, double *states, const double *parameters,
                       double *algebraicsForTransfer, const int *algebraicsForTransferIndices, int nAlgebraicsForTransferIndices,
                       const double *elementLengths, char *firingEvents, int firingEventsNRows, int firingEventsNColumns,
                       double *setSpecificStatesFrequencyJitter, int frequencyJitterNColumns, char *fiberIsCurrentlyStimulated,
                       double startTime, double timeStepWidthSplitting, int nTimeStepsSplitting, double dt0D, int nTimeSteps0D, double dt1D, int nTimeSteps1D,
                       double prefactor, double valueForStimulatedPoint)
{
)";

  if (optimizationType_ == "gpu")
    sourceCode << R"(
  #pragma omp target map(tofrom: fiberData, states) map(to: parameters, algebraicsForTransferIndices, elementLengths, firingEvents, setSpecificStatesFrequencyJitter) map(from: algebraicsForTransfer)
  {
)";
   sourceCode << R"(
  // size constants
  const int nInstancesPerFiber = )" << nInstancesToComputePerFiber_ << R"(;
  const int nElementsOnFiber = )" << nInstancesToComputePerFiber_-1 << R"(;
  const int nFibersToCompute = )" << nFibersToCompute_ << R"(;
  const long long nInstancesToCompute = )" << nFibersToCompute_*nInstancesToComputePerFiber_ << R"(;  // = nInstancesPerFiber*nFibersToCompute
  //const int nStates = )" << nStates << R"(;
  //const int nAlgebraics = )" << nAlgebraics << R"(;

  // loop over splitting time steps
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
    sourceCode << "\n    #pragma omp teams distribute parallel for collapse(2)";
  sourceCode << R"(
    for (int fiberNo = 0; fiberNo < nFibersToCompute; fiberNo++)
    {
      // loop over instances to compute here
      for (int instanceNo = 0; instanceNo < nInstancesPerFiber; instanceNo++)
      {
        int instanceToComputeNo = fiberNo*nInstancesPerFiber + instanceNo;    // index of instance over all fibers

        // determine if current point is at center of fiber
        int fiberCenterIndex = fiberData[fiberNo].fiberStimulationPointIndex;
        bool currentPointIsInCenter = fabs(fiberCenterIndex - instanceNo) < 4;

        // loop over 0D timesteps
        for (int timeStepNo = 0; timeStepNo < nTimeSteps0D; timeStepNo++)
        {
          double currentTime = currentTimeSplitting + timeStepNo * dt0D;

          // determine if fiber gets stimulated
          // check if current point will be stimulated
          bool stimulateCurrentPoint = false;
          if (currentPointIsInCenter)
            stimulateCurrentPoint = isCurrentPointStimulated(
              fiberData, firingEvents, firingEventsNRows, firingEventsNColumns,
              setSpecificStatesFrequencyJitter, frequencyJitterNColumns, fiberIsCurrentlyStimulated,
              fiberNo, currentTime, currentPointIsInCenter);
          const bool storeAlgebraicsForTransfer = false;

          // call method to compute 0D problem
          // compute0DInstance_(fiberPointBuffers_[pointBuffersNo].states, fiberPointBuffersParameters_[pointBuffersNo],
          //                  currentTime, timeStepWidth, stimulateCurrentPoint,
          //                  argumentStoreAlgebraics, fiberPointBuffersAlgebraicsForTransfer_[pointBuffersNo],
          //                  algebraicsForTransferIndices_, valueForStimulatedPoint_);
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
    #pragma omp teams distribute parallel for
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

        u_center = states[0 + instanceToComputeNo];  // state 0 of the current instance

        // contribution from left element
        if (valueNo > 0)
        {
          u_previous = states[0 + instanceToComputeNo - 1];  // state 0 of the left instance

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
          u_next = states[0 + instanceToComputeNo + 1];  // state 0 of the right instance

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
      states[0+nValues-1] = dIntermediate[nValues-1];  // state 0 of the point (nValues-1)

      double previousValue = dIntermediate[nValues-1];

      // loop over entries / rows of matrices
      for (int valueNo = nValues-2; valueNo >= 0; valueNo--)
      {
        int instanceToComputeNo = fiberNo*nInstancesPerFiber + valueNo;

        // x_i = d'_i - c'_i * x_{i+1}
        double resultValue = dIntermediate[valueNo] - cIntermediate[valueNo] * previousValue;
        states[0+instanceToComputeNo] = resultValue;

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

    #pragma omp teams distribute parallel for collapse(2)
    for (int fiberNo = 0; fiberNo < nFibersToCompute; fiberNo++)
    {)";
  sourceCode << R"(
      // loop over instances to compute here
      for (int instanceNo = 0; instanceNo < nInstancesPerFiber; instanceNo++)
      {
        int instanceToComputeNo = fiberNo*nInstancesPerFiber + instanceNo;    // index of instance over all fibers

        // determine if current point is at center of fiber
        int fiberCenterIndex = fiberData[fiberNo].fiberStimulationPointIndex;
        bool currentPointIsInCenter = fabs(fiberCenterIndex - instanceNo) < 4;

        // loop over 0D timesteps
        for (int timeStepNo = 0; timeStepNo < nTimeSteps0D; timeStepNo++)
        {
          double currentTime = midTimeSplitting + timeStepNo * dt0D;

          // determine if fiber gets stimulated
          // check if current point will be stimulated
          bool stimulateCurrentPoint = false;
          if (currentPointIsInCenter)
            stimulateCurrentPoint = isCurrentPointStimulated(
              fiberData, firingEvents, firingEventsNRows, firingEventsNColumns,
              setSpecificStatesFrequencyJitter, frequencyJitterNColumns, fiberIsCurrentlyStimulated,
              fiberNo, currentTime, currentPointIsInCenter);
          const bool storeAlgebraicsForTransfer = storeAlgebraicsForTransferSplitting && timeStepNo == nTimeSteps0D-1;

          // call method to compute 0D problem
          // compute0DInstance_(fiberPointBuffers_[pointBuffersNo].states, fiberPointBuffersParameters_[pointBuffersNo],
          //                  currentTime, timeStepWidth, stimulateCurrentPoint,
          //                  argumentStoreAlgebraics, fiberPointBuffersAlgebraicsForTransfer_[pointBuffersNo],
          //                  algebraicsForTransferIndices_, valueForStimulatedPoint_);

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

  const long long int nInstancesToCompute = nFibersToCompute_*nInstancesToComputePerFiber_;

  // copy all states values from fiberData to the gpuStates_ vector which is then used in the generated code
  for (int fiberDataNo = 0; fiberDataNo < nFibersToCompute_; fiberDataNo++)
  {
    for (int instanceNo = 0; instanceNo < nInstancesToComputePerFiber_; instanceNo++)
    {
      const int instanceNoTotal = fiberDataNo*nInstancesToComputePerFiber_ + instanceNo;

      // set state 0
      gpuStates_[0*nInstancesToCompute + instanceNoTotal] = fiberData_[fiberDataNo].vmValues[instanceNo];
    }
  }

  // call the compiled function
  computeMonodomain_(fiberData_.data(), gpuStates_.data(), gpuParameters_.data(),
                     gpuAlgebraicsForTransfer_.data(), algebraicsForTransferIndices_.data(), algebraicsForTransferIndices_.size(),
                     gpuElementLengths_.data(), gpuFiringEvents_.data(), gpuFiringEventsNRows_, gpuFiringEventsNColumns_,
                     gpuSetSpecificStatesFrequencyJitter_.data(), gpuFrequencyJitterNColumns_, gpuFiberIsCurrentlyStimulated_.data(),
                     startTime, timeStepWidthSplitting, nTimeStepsSplitting_, dt0D, nTimeSteps0D, dt1D, nTimeSteps1D,
                     prefactor, valueForStimulatedPoint_);

  // copy the resulting values back to fiberData_
  for (int fiberDataNo = 0; fiberDataNo < nFibersToCompute_; fiberDataNo++)
  {
    for (int instanceNo = 0; instanceNo < nInstancesToComputePerFiber_; instanceNo++)
    {
      const int instanceNoTotal = fiberDataNo*nInstancesToComputePerFiber_ + instanceNo;

      assert(statesForTransferIndices_.size() > 0);
      const int stateToTransfer = statesForTransferIndices_[0];  // transfer the first state value
      fiberData_[fiberDataNo].vmValues[instanceNo] = gpuStates_[stateToTransfer*nInstancesToCompute + instanceNoTotal];

      // loop over further states to transfer
      int furtherDataIndex = 0;
      for (int i = 1; i < statesForTransferIndices_.size(); i++, furtherDataIndex++)
      {
        const int stateToTransfer = statesForTransferIndices_[i];

        fiberData_[fiberDataNo].furtherStatesAndAlgebraicsValues[furtherDataIndex*nInstancesToComputePerFiber_ + instanceNo]
          = gpuStates_[stateToTransfer*nInstancesToCompute + instanceNoTotal];
      }

      // loop over algebraics to transfer
      for (int i = 0; i < algebraicsForTransferIndices_.size(); i++, furtherDataIndex++)
      {
        fiberData_[fiberDataNo].furtherStatesAndAlgebraicsValues[furtherDataIndex*nInstancesToComputePerFiber_ + instanceNo]
          = gpuAlgebraicsForTransfer_[i*nInstancesToCompute + instanceNoTotal];
      }
    }
    LOG(DEBUG) << "states and algebraics for transfer at fiberDataNo=" << fiberDataNo << ": " << fiberData_[fiberDataNo].furtherStatesAndAlgebraicsValues;
    LOG(DEBUG) << "size: " << fiberData_[fiberDataNo].furtherStatesAndAlgebraicsValues.size() << ", nInstancesToComputePerFiber_: " << nInstancesToComputePerFiber_;
  }


  currentTime_ = instances[0].endTime();
}
