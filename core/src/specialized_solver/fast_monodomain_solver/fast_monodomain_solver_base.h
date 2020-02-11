#pragma once

#include <Python.h>  // has to be the first included header
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

/** The implementation of a monodomain solver as used in the fibers_emg example, number of states and intermediates is templated.
 *  This class contains all functionality except the reaction term. Deriving classes only need to implement compute0D.
  */
template<int nStates, int nIntermediates>
class FastMonodomainSolverBase : public Runnable
{
public:

  typedef FunctionSpace::FunctionSpace<
    Mesh::StructuredDeformableOfDimension<1>,
    BasisFunction::LagrangeOfOrder<1>
  > FiberFunctionSpace;

  typedef CellmlAdapter<
    nStates, nIntermediates,  // nStates,nIntermediates: 57,1 = Shorten, 4,9 = Hodgkin Huxley
    FiberFunctionSpace
  > CellmlAdapterType;

  typedef TimeSteppingScheme::ImplicitEuler<          // fiber diffusion, note that implicit euler gives lower error in this case than crank nicolson
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredDeformableOfDimension<1>,
      BasisFunction::LagrangeOfOrder<1>,
      Quadrature::Gauss<2>,
      Equation::Dynamic::IsotropicDiffusion
    >
  > ImplicitEuler;

  typedef Control::MultipleInstances<                       // fibers
    OperatorSplitting::Strang<
      Control::MultipleInstances<
        TimeSteppingScheme::Heun<                   // fiber reaction term
          CellmlAdapterType
        >
      >,
      Control::MultipleInstances<
        ImplicitEuler
      >
    >
  > NestedSolversType;

  typedef typename NestedSolversType::FunctionSpace FunctionSpace;
  typedef typename NestedSolversType::Data Data;
  typedef typename NestedSolversType::OutputConnectorDataType OutputConnectorDataType;

  //! constructor
  FastMonodomainSolverBase(const DihuContext &context);

  //! initialize the simulation, this is called from run
  void initialize();

  //! run the whole simulation
  void run();

  //! run simulation for the specified timespan
  void advanceTimeSpan();

  //! reset to uninitialized state
  void reset();

  //! get a reference to the data object
  Data &data();

  //! set a new time interval that will be simulated by next call to advanceTimeSpan.
  void setTimeSpan(double startTime, double endTime);

  //! get the output connector data, to be used for a surrounding solver
  std::shared_ptr<OutputConnectorDataType> getOutputConnectorData();

protected:

  //! create a source file with compute0D function from the CellML model
  void initializeCellMLSourceFile();

  //! get element lengths and vmValues from the other ranks
  void fetchFiberData();

  //! send vmValues data from fiberData_ back to the fibers where it belongs to and set in the respective field variable
  void updateFiberData();

  //! solve the 0D problem, starting from startTime. This is the part that is usually provided by the cellml file
  void compute0D(double startTime, double timeStepWidth, int nTimeSteps, bool storeIntermediatesForTransfer);

  //! compute one time step of the right hand side for a single simd vector of instances
  virtual void compute0DInstance(Vc::double_v states[], std::vector<Vc::double_v> &parameters, double currentTime, double timeStepWidth,
                                 bool stimulate, bool storeIntermediatesForTransfer,
                                 std::vector<Vc::double_v> &intermediatesForTransfer){};

  //! solve the 1D problem (diffusion), starting from startTime
  void compute1D(double startTime, double timeStepWidth, int nTimeSteps, double prefactor);

  //! compute the 0D-1D problem with Strang splitting
  void computeMonodomain();

  //! set the initial values for all states
  virtual void initializeStates(Vc::double_v states[]){};

  PythonConfig specificSettings_;    //< config for this object

  NestedSolversType nestedSolvers_;   //< the nested solvers object that would normally solve the problem

  /** data to be exchanged for computation of a single fiber
   *  The data stored herein is used for local computation.
   */
  struct FiberData
  {
    std::vector<double> elementLengths;   //< lengths of the 1D elements
    std::vector<double> vmValues;         //< values of Vm
    std::vector<double> furtherStatesAndIntermediatesValues;    //< all data to be transferred back to the fibers, apart from vmValues, corresponding to statesForTransfer_ and intermediatesForTransfer_ (array of struct memory layout)
    int valuesLength;                     //< number of vmValues
    global_no_t valuesOffset;             //< number of vmValues in previous entries in fiberData_

    int fiberNoGlobal;                    //< fiberNo as given in settings (value of additionalArgument)
    int motorUnitNo;                      //< motor unit no.

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

  OutputWriter::Manager outputWriterManager_;     ///< manager object holding all output writers

  std::vector<FiberData> fiberData_;  //< vector of fibers,
  int nFibersToCompute_;              //< number of fibers where own rank is involved (>= n.fibers that are computed by own rank)
  int nInstancesToCompute_;           //< number of instances of the Hodgkin-Huxley problem to compute on this rank
  double currentTime_;                //< the current time used for the output writer
  int nTimeStepsSplitting_;           //< number of times to repeat the Strang splitting for one advanceTimeSpan() call of FastMonodomainSolver

  std::vector<int> statesForTransfer_;                          //< state no.s to transfer to other solvers within output connector data
  std::vector<int> intermediatesForTransfer_;   //< which intermediates should be transferred to other solvers as part of output connector data

  std::vector<std::vector<Vc::double_v>> fiberPointBuffersParameters_;        //< constant parameter values, changing parameters is not implemented
  std::vector<std::vector<Vc::double_v>> fiberPointBuffersIntermediatesForTransfer_;   //<  [fiberPointNo][intermediateToTransferNo], intermediate values to use for output connector data

  CellmlSourceCodeGenerator cellmlSourceCodeGenerator_;    ///< object that holds all source code related to the model

  void (*compute0DInstance_)(Vc::double_v [], std::vector<Vc::double_v> &, double, double, bool, bool, std::vector<Vc::double_v> &, const std::vector<int> &);   //< runtime-created and loaded function to compute one Heun step of the 0D problem
  void (*initializeStates_)(Vc::double_v states[]);  //< runtime-created and loaded function to set all initial values for the states

};

#include "specialized_solver/fast_monodomain_solver/fast_monodomain_solver_base.tpp"
