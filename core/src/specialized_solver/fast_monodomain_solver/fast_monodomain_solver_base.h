#pragma once

#include <Python.h>  // has to be the first included header
#include <array>
#include <Vc/Vc>

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

  //! get the slot connector data, to be used for an enclosing solver
  std::shared_ptr<SlotConnectorDataType> getSlotConnectorData();

  //! get a reference to the nested solvers
  NestedSolversType &nestedSolvers();

protected:

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
  struct FiberData;
  bool isCurrentPointStimulated(int fiberDataNo, double currentTime, bool currentPointIsInCenter);

  //! method to be called after the compute0D, updates the information in fiberPointBuffersStatesAreCloseToEquilibrium_
  void equilibriumAccelerationUpdate(const Vc::double_v statesPreviousValues[], int pointBuffersNo);

  //! check if the 0D computations for the current point are disabled because the states are in equilibrium
  bool isEquilibriumAccelerationCurrentPointDisabled(bool stimulateCurrentPoint, int pointBuffersNo);

  //! set the initial values for all states
  virtual void initializeStates(Vc::double_v states[]){};

  //! generate source code to solve the monodomain equation with offloading to gpu
  void addMonodomainSolverGpuSource(std::string filename);

  //! call the GPU code to compute the monodomain equation for advanceTimestep
  void computeMonodomainGpu();

  PythonConfig specificSettings_;    //< config for this object

  NestedSolversType nestedSolvers_;   //< the nested solvers object that would normally solve the problem

  /** data to be exchanged for computation of a single fiber
   *  The data stored herein is used for local computation.
   */
  struct FiberData
  {
    std::vector<double> elementLengths;   //< lengths of the 1D elements
    std::vector<double> vmValues;         //< values of Vm
    std::vector<double> furtherStatesAndAlgebraicsValues;    //< all data to be transferred back to the fibers, apart from vmValues, corresponding to statesForTransfer_ and algebraicsForTransfer_ (array of struct memory layout)
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

  std::vector<FiberPointBuffers<nStates>> fiberPointBuffers_;    //< computation buffers for the 0D problem

  std::string fiberDistributionFilename_;  //< filename of the fiberDistributionFile, which contains motor unit numbers for fiber numbers
  std::string firingTimesFilename_;        //< filename of the firingTimesFile, which contains points in time of stimulation for each motor unit

  std::vector<std::vector<bool>> firingEvents_;   //< if a motor unit firingEvents_[timeStepNo][motorUnitNo]
  std::vector<int> motorUnitNo_;                  //< number of motor unit for given fiber no motorUnitNo_[fiberNo]
  std::string durationLogKey0D_;                  //< duration log key for the 0D problem
  std::string durationLogKey1D_;                  //< duration log key for the 1D problem

  OutputWriter::Manager outputWriterManager_;     //< manager object holding all output writers

  std::vector<FiberData> fiberData_;  //< vector of fibers, the number of entries is the number of fibers to be computed by the own rank (nFibersToCompute_)
  int nFibersToCompute_;              //< number of fibers where own rank is involved (>= n.fibers that are computed by own rank)
  int nInstancesToCompute_;           //< number of instances of the Hodgkin-Huxley problem to compute on this rank
  double currentTime_;                //< the current time used for the output writer
  int nTimeStepsSplitting_;           //< number of times to repeat the Strang splitting for one advanceTimeSpan() call of FastMonodomainSolver

  bool onlyComputeIfHasBeenStimulated_;       //< option if fiber should only be computed after it has been stimulated for the first time
  std::vector<bool> fiberHasBeenStimulated_;  //< for every fiber if it has been stimulated

  bool disableComputationWhenStatesAreCloseToEquilibrium_;                  //< option to avoid computation when the states won't change much
  enum state_t {
    constant,                         //< the state values at the own point did not change in the last computation (according to a tolerance). This means the current point does not need to be computed.
    neighbour_not_constant,           //< the state values did not change, so the state is constant, but at a neighbouring point the value changed. This means the own value has to be computed because it can change due to diffusion.
    not_constant                      //< the state values at the own point change and have to be computed
  };                                                                        //< type for fiberPointBuffersStatesAreCloseToEquilibrium_
  std::vector<state_t> fiberPointBuffersStatesAreCloseToEquilibrium_;       //< for every entry in fiberPointBuffers_, constant if the states didn't change too much in the last compute0D, neighbour_not_constant if the state of the neighbouring pointBuffer changes
  int nFiberPointBufferStatesCloseToEquilibrium_;                           //< number of "constant" entries in fiberPointBuffersStatesAreCloseToEquilibrium_

  std::vector<int> statesForTransfer_;          //< state no.s to transfer to other solvers within slot connector data
  std::vector<int> algebraicsForTransfer_;      //< which algebraics should be transferred to other solvers as part of slot connector data
  std::vector<double> parameters_;              //< parameters vector
  double valueForStimulatedPoint_;              //< value to which the first state will be set if stimulated
  double neuromuscularJunctionRelativeSize_;    //< relative size of the range where the neuromuscular junction is located

  std::vector<std::vector<Vc::double_v>> fiberPointBuffersParameters_;        //< constant parameter values, changing parameters is not implemented
  std::vector<std::vector<Vc::double_v>> fiberPointBuffersAlgebraicsForTransfer_;   //<  [fiberPointNo][algebraicToTransferNo], algebraic values to use for slot connector data

  void (*compute0DInstance_)(Vc::double_v [], std::vector<Vc::double_v> &, double, double, bool, bool, std::vector<Vc::double_v> &, const std::vector<int> &, double);   //< runtime-created and loaded function to compute one Heun step of the 0D problem
  void (*initializeStates_)(Vc::double_v states[]);  //< runtime-created and loaded function to set all initial values for the states

  bool useVc_;                                       //< if the Vc library is used, if not, code for the GPU is generated
  bool initialized_;                                 //< if initialize was already called
};

#include "specialized_solver/fast_monodomain_solver/fast_monodomain_solver_base.tpp"
#include "specialized_solver/fast_monodomain_solver/fast_monodomain_solver_communication.tpp"
#include "specialized_solver/fast_monodomain_solver/fast_monodomain_solver_compute.tpp"
#include "specialized_solver/fast_monodomain_solver/fast_monodomain_solver_initialization.tpp"
