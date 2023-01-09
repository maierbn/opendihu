#pragma once

#include <Python.h>  // has to be the first included header
#include <array>
#include <vc_or_std_simd.h>  // this includes <Vc/Vc> or a Vc-emulating wrapper of <experimental/simd> if available

#include "control/multiple_instances.h"
#include "operator_splitting/strang.h"
#include "time_stepping_scheme/heun.h"
#include "cellml/03_cellml_adapter.h"
#include "function_space/function_space.h"
#include "mesh/structured_deformable.h"
#include "basis_function/lagrange.h"
#include "time_stepping_scheme/implicit_euler.h"
#include "spatial_discretization/finite_element_method/finite_element_method.h"

/** Buffers for CellML computation
  *  Includes Vc::double_v::size() instances of the CellML problem (usually 4 when using AVX-2).
  *  These will be computed at once using vector instructions.
  */
template<int nStates>
struct FiberPointBuffers
{
  Vc::double_v states[nStates];
};

#ifndef HAVE_STDSIMD
// only if we are using Vc, it is not necessary for std::simd

/** Specialize the default allocator for the FiberPointBuffers struct to use the aligned allocated provided by Vc.
 *  This could also be done by Vc_DECLARE_ALLOCATOR(<class>), but not here because of the template parameter nStates.
 */
namespace std
{
template<int nStates>
class allocator<FiberPointBuffers<nStates>> :
  public ::Vc::Allocator<FiberPointBuffers<nStates>>
{
public:
  template <typename U>
  struct rebind
  {
    typedef ::std::allocator<U> other;
  };
};
}
#endif

/** The implementation of a monodomain solver as used in the fibers_emg example, number of states and algebraics is templated.
 *  This class contains all functionality except the reaction term. Deriving classes only need to implement compute0D.
  */
template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
class FastMonodomainSolverBase : public Runnable
{
public:

  typedef FunctionSpace::FunctionSpace<
    Mesh::StructuredDeformableOfDimension<1>,
    BasisFunction::LagrangeOfOrder<1>
  > FiberFunctionSpace;

  typedef CellmlAdapter<
    nStates, nAlgebraics,  // nStates,nAlgebraics: 57,1 = Shorten, 4,9 = Hodgkin Huxley
    FiberFunctionSpace
  > CellmlAdapterType;

  typedef Control::MultipleInstances<                       // fibers
    OperatorSplitting::Strang<
      Control::MultipleInstances<
        TimeSteppingScheme::Heun<                   // fiber reaction term
          CellmlAdapterType
        >
      >,
      Control::MultipleInstances<
        DiffusionTimeSteppingScheme
      >
    >
  > NestedSolversType;

  typedef typename NestedSolversType::FunctionSpace FunctionSpace;
  typedef typename NestedSolversType::Data Data;
  typedef typename NestedSolversType::SlotConnectorDataType SlotConnectorDataType;

  //! constructor
  FastMonodomainSolverBase(const DihuContext &context);

  //! initialize the simulation, this is called from run
  void initialize();

  //! run the whole simulation
  void run();

  //! run simulation for the specified timespan
  void advanceTimeSpan(bool withOutputWritersEnabled = true);

  //! call the output writer on the data object, output files will contain currentTime, with callCountIncrement !=1 output timesteps can be skipped
  void callOutputWriter(int timeStepNo, double currentTime, int callCountIncrement = 1);

  //! reset to uninitialized state
  void reset();

  //! get a reference to the data object
  Data &data();

  //! set a new time interval that will be simulated by next call to advanceTimeSpan.
  void setTimeSpan(double startTime, double endTime);

  //! start time of time interval to be simulated
  double startTime();

  //! end time of simulation
  double endTime();

  //! get the slot connector data, to be used for an enclosing solver
  std::shared_ptr<SlotConnectorDataType> getSlotConnectorData();

  //! get a reference to the nested solvers
  NestedSolversType &nestedSolvers();

  /** data to be exchanged for computation of a single fiber
   *  The data stored herein is used for local computation.
   */
  struct FiberData
  {
    std::vector<double> elementLengths;   //< lengths of the 1D elements
    std::vector<double> vmValues;         //< values of Vm
    std::vector<double> furtherStatesAndAlgebraicsValues;    //< all data to be transferred back to the fibers, apart from vmValues, corresponding to statesForTransferIndices_ and algebraicsForTransferIndices_ (array of struct memory layout)
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

  void updateFiberState();

  void saveFiberState();



protected:

  //! load the firing times file and initialize the firingEvents_ and motorUnitNo_ variables
  void initializeFiringTimes();

  //! initialize the internal data structures, such as fiberData_
  void initializeDataStructures();

  //! set the names of the field variables in the data connector slots
  void initializeFieldVariableNames();

  //! create a source file with compute0D function from the CellML model, using the vc optimization type
  void initializeCellMLSourceFileVc();

  //! create a source file with compute0D function from the CellML model, using the gpu optimization type
  void initializeCellMLSourceFileGpu();

  //! get element lengths and vmValues from the other ranks
  void fetchFiberData();

  //! send vmValues data from fiberData_ back to the fibers where it belongs to and set in the respective field variable
  void updateFiberData();

  //! solve the 0D problem, starting from startTime. This is the part that is usually provided by the cellml file
  void compute0D(double startTime, double timeStepWidth, int nTimeSteps, bool storeAlgebraicsForTransfer);

  //! compute one time step of the right hand side for a single simd vector of instances
  virtual void compute0DInstance(Vc::double_v states[], std::vector<Vc::double_v> &parameters, double currentTime, double timeStepWidth,
                                 bool stimulate, bool storeAlgebraicsForTransfer,
                                 std::vector<Vc::double_v> &algebraicsForTransfer){};

  //! solve the 1D problem (diffusion), starting from startTime
  void compute1D(double startTime, double timeStepWidth, int nTimeSteps, double prefactor);

  //! compute the 0D-1D problem with Strang splitting
  void computeMonodomain();

  //! check if the current point will be stimulated now
  bool isCurrentPointStimulated(int fiberDataNo, double currentTime, bool currentPointIsInCenter);

  //! method to be called after the compute0D, updates the information in fiberPointBuffersStatesAreCloseToEquilibrium_
  void equilibriumAccelerationUpdate(const Vc::double_v statesPreviousValues[], int pointBuffersNo);

  //! check if the 0D computations for the current point are disabled because the states are in equilibrium
  bool isEquilibriumAccelerationCurrentPointDisabled(bool stimulateCurrentPoint, int pointBuffersNo);

  //! set the initial values for all states
  virtual void initializeStates(Vc::double_v states[]){};

  //! initialize the states vector and other static data that is used for GPU computation
  void initializeValuesOnGpu();

  //! generate source code to solve the monodomain equation with offloading to gpu
  void generateMonodomainSolverGpuSource(std::string filename, std::string headerCode, std::string mainCode);

  //! call the GPU code to compute the monodomain equation for advanceTimestep
  void computeMonodomainGpu();
  
  //! find out the currently used version of GCC, return as a string in the format, e.g., "10.2.0"
  std::string checkGccVersion();

  PythonConfig specificSettings_;    //< config for this object

  NestedSolversType nestedSolvers_;   //< the nested solvers object that would normally solve the problem

  std::vector<FiberPointBuffers<nStates>> fiberPointBuffers_;    //< computation buffers for the 0D problem, the states vector used when optimizationType == "vc"

  std::string fiberDistributionFilename_;  //< filename of the fiberDistributionFile, which contains motor unit numbers for fiber numbers
  std::string firingTimesFilename_;        //< filename of the firingTimesFile, which contains points in time of stimulation for each motor unit

  std::vector<std::vector<bool>> firingEvents_;   //< if a motor unit fires, firingEvents_[timeStepNo][motorUnitNo]
  std::vector<int> motorUnitNo_;                  //< number of motor unit for given fiber no motorUnitNo_[fiberNo]
  std::string durationLogKey0D_;                  //< duration log key for the 0D problem
  std::string durationLogKey1D_;                  //< duration log key for the 1D problem

  OutputWriter::Manager outputWriterManager_;     //< manager object holding all output writers
  std::vector<FiberData> fiberDataOld_;
  std::vector<FiberData> fiberData_;  //< vector of fibers, the number of entries is the number of fibers to be computed by the own rank (nFibersToCompute_)
  int nFibersToCompute_;              //< number of fibers where own rank is involved (>= n.fibers that are computed by own rank)
  int nInstancesToCompute_;           //< number of instances of the Hodgkin-Huxley (or other CellML) problem to compute on this rank
  int nInstancesToComputePerFiber_;   //< number of instances to compute per fiber, i.e., global number of instances of a fiber
  int nParametersPerInstance_;        //< number of parameters per instance
  double currentTime_;                //< the current time used for the output writer
  int nTimeStepsSplitting_;           //< number of times to repeat the Strang splitting for one advanceTimeSpan() call of FastMonodomainSolver

  bool onlyComputeIfHasBeenStimulated_;       //< option if fiber should only be computed after it has been stimulated for the first time
  std::vector<bool> fiberHasBeenStimulated_;  //< for every fiber if it has been stimulated

  bool disableComputationWhenStatesAreCloseToEquilibrium_;                  //< option to avoid computation when the states won't change much
  enum state_t {
    inactive,                         //< the state values at the own point did not change in the last computation (according to a tolerance). This means the current point does not need to be computed.
    neighbor_is_active,           //< the state values did not change, so the state is inactive, but at a neighbouring point the value changed. This means the own value has to be computed because it can change due to diffusion.
    active                      //< the state values at the own point change and have to be computed
  };                                                                        //< type for fiberPointBuffersStatesAreCloseToEquilibrium_
  std::vector<state_t> fiberPointBuffersStatesAreCloseToEquilibrium_;       //< for every entry in fiberPointBuffers_, inactive if the states didn't change too much in the last compute0D, neighbor_is_active if the state of the neighbouring pointBuffer changes
  int nFiberPointBufferStatesCloseToEquilibrium_;                           //< number of "inactive" entries in fiberPointBuffersStatesAreCloseToEquilibrium_
  bool setComputeStateInformation_;                                         //< whether the information in fiberPointBuffersStatesAreCloseToEquilibrium_ should be added to the algebraics to transfer in a variable named "computeStateInformation"

  std::vector<int> statesForTransferIndices_;          //< state no.s to transfer to other solvers within slot connector data
  std::vector<int> algebraicsForTransferIndices_;      //< which algebraics should be transferred to other solvers as part of slot connector data
  double valueForStimulatedPoint_;              //< value to which the first state will be set if stimulated
  double neuromuscularJunctionRelativeSize_;    //< relative size of the range where the neuromuscular junction is located

  std::vector<std::vector<Vc::double_v>> fiberPointBuffersParameters_;        //< constant parameter values, changing parameters is not implemented
  std::vector<std::vector<Vc::double_v>> fiberPointBuffersAlgebraicsForTransfer_;   //<  [fiberPointNo][algebraicToTransferNo], algebraic values to use for slot connector data

  std::vector<float> gpuParameters_;              //< for "gpu": constant parameter values, in struct of array memory layout: gpuParameters_[parameterNo*nInstances + instanceNo]
  std::vector<double> gpuAlgebraicsForTransfer_;   //< for "gpu": algebraic values to use for slot connector data, in struct of array memory layout: gpuAlgebraicsForTransfer_[algebraicNo*nInstances + instanceNo]
  std::vector<double> gpuStatesForTransfer_;       //< for "gpu": state values to use for slot connector data, in struct of array memory layout: gpuStatesForTransfer_[stateInThisListIndex*nInstances + instanceNo]
  std::vector<float> gpuElementLengths_;          //< for "gpu": the lengths of the 1D elements, in struct of array memory layout: gpuElementLengths_[fiberDataNo*nElementsOnFiber + elementNo]
  std::vector<char> gpuFiringEvents_;              //< for "gpu": if a motor unit fires at a specified time, 1=yes, 0=no, gpuFiringEvents_[timeStepNo*nMotorUnits + motorUnitNo]
  std::vector<double> gpuSetSpecificStatesFrequencyJitter_;  //< for "gpu", value of option with the same name in the python settings: gpuSetSpecificStatesFrequencyJitter_[fiberNo*nColumns + columnNo]
  int gpuFiringEventsNRows_;                       //< for "gpu": number of rows in the firing events file
  int gpuFiringEventsNColumns_;                    //< for "gpu": number of columns in the firing events file
  int gpuFrequencyJitterNColumns_;                 //< for "gpu": number of columns in the gpuSetSpecificStatesFrequencyJitter_ array
  std::vector<char> gpuFiberIsCurrentlyStimulated_; //< for "gpu": the value of fiberData_[].currentlyStimulating
  std::vector<int> gpuMotorUnitNo_;                      //< motor unit no.
  std::vector<int> gpuFiberStimulationPointIndex_;       //< index of the point on the fiber where to stimulate, i.e. position of the neuromuscular junction, if at center, it is equal to (int)(fiberData_[fiberDataNo].valuesLength / 2)
  std::vector<double> gpuLastStimulationCheckTime_;      //< last time the fiber was checked for stimulation
  std::vector<double> gpuSetSpecificStatesCallFrequency_;        //< value of option with the same name in the python settings
  std::vector<double> gpuSetSpecificStatesRepeatAfterFirstCall_; //< how long in ms the prescribed value should be set
  std::vector<double> gpuSetSpecificStatesCallEnableBegin_;      //< value of option with the same name in the python settings
  std::vector<double> gpuCurrentJitter_;                         //< current absolute value of jitter to add to setSpecificStatesCallFrequency
  std::vector<int> gpuJitterIndex_;                              //< index of the vector in setSpecificStatesFrequencyJitter which is the current value to use
  std::vector<double> gpuVmValues_;                      //< for "gpu": values of the first state, gpuVmValues_[instanceToComputeNo]
  bool generateGpuSource_;                               //< if the GPU source code should be generated, if not it reuses the existing file, this is for debugging

  void (*compute0DInstance_)(Vc::double_v [], std::vector<Vc::double_v> &, double, double, bool, bool, std::vector<Vc::double_v> &, const std::vector<int> &, double);   //< runtime-created and loaded function to compute one Heun step of the 0D problem
  void (*computeMonodomain_)(const float *parameters,
                              double *algebraicsForTransfer, double *statesForTransfer, const float *elementLengths,
                              double startTime, double timeStepWidthSplitting, int nTimeStepsSplitting, double dt0D, int nTimeSteps0D, double dt1D, int nTimeSteps1D,
                              double prefactor, double valueForStimulatedPoint);   //< runtime-created and loaded function to compute monodomain equation
  void (*initializeArrays_)(const double *statesOneInstance, const int *algebraicsForTransferIndicesParameter, const int *statesForTransferIndicesParameter,
                            const char *firingEventsParameter, const double *setSpecificStatesFrequencyJitterParameter, const int *motorUnitNoParameter,
                            const int *fiberStimulationPointIndexParameter, const double *lastStimulationCheckTimeParameter,
                            const double *setSpecificStatesCallFrequencyParameter, const double *setSpecificStatesRepeatAfterFirstCallParameter,
                            const double *setSpecificStatesCallEnableBeginParameter);   //< function that initializes all data on the target device (GPU)
  void (*initializeStates_)(Vc::double_v states[]);  //< runtime-created and loaded function to set all initial values for the states

  bool useVc_;                                       //< if the Vc library is used, if not, code for the GPU or OpenMP is generated
  std::string optimizationType_;                     //< the optimization type as given in the settings, one of "vc", "openmp", "gpu"
  bool initialized_;                                 //< if initialize was already called
};

#include "specialized_solver/fast_monodomain_solver/fast_monodomain_solver_base.tpp"
#include "specialized_solver/fast_monodomain_solver/fast_monodomain_solver_communication.tpp"
#include "specialized_solver/fast_monodomain_solver/fast_monodomain_solver_compute.tpp"
#include "specialized_solver/fast_monodomain_solver/fast_monodomain_solver_initialization.tpp"
#include "specialized_solver/fast_monodomain_solver/fast_monodomain_solver_gpu.tpp"
