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

/** Buffers for Hodgkin-Huxley computation
  *  Includes Vc::double_v::size() instances of the HH problem (usually 4).
  *  These will be computed at once using vector instructions.
  */
struct FiberPointBuffers
{
  Vc::double_v states[4];
};
Vc_DECLARE_ALLOCATOR(FiberPointBuffers)

/** A specialized implementation for the monodomain equation, as in fibers_emg example.
 *  This is the general class, relevant is the specialization below.
  */
template<typename FibersEMG>
class FastMonodomainSolver
{
};

/** Partial specialization for exactly the combination of nested solvers for which it is optimized.
 *  (fibers_emg for Hodgkin-Huxley model)
 */
template<>
class FastMonodomainSolver<
  Control::MultipleInstances<                       // fibers
    OperatorSplitting::Strang<
      Control::MultipleInstances<
        TimeSteppingScheme::Heun<                   // fiber reaction term
          CellmlAdapter<
            4,9,  // nStates,nIntermediates: 57,1 = Shorten, 4,9 = Hodgkin Huxley
            FunctionSpace::FunctionSpace<
              Mesh::StructuredDeformableOfDimension<1>,
              BasisFunction::LagrangeOfOrder<1>
            >
          >
        >
      >,
      Control::MultipleInstances<
        TimeSteppingScheme::ImplicitEuler<          // fiber diffusion, note that implicit euler gives lower error in this case than crank nicolson
          SpatialDiscretization::FiniteElementMethod<
            Mesh::StructuredDeformableOfDimension<1>,
            BasisFunction::LagrangeOfOrder<1>,
            Quadrature::Gauss<2>,
            Equation::Dynamic::IsotropicDiffusion
          >
        >
      >
    >
  >
> : public Runnable
{
public:

  typedef FunctionSpace::FunctionSpace<
    Mesh::StructuredDeformableOfDimension<1>,
    BasisFunction::LagrangeOfOrder<1>
  > FiberFunctionSpace;

  typedef CellmlAdapter<
    4,9,  // nStates,nIntermediates: 57,1 = Shorten, 4,9 = Hodgkin Huxley
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
  FastMonodomainSolver(const DihuContext &context);

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
  OutputConnectorDataType getOutputConnectorData();

private:

  //! get element lengths and vmValues from the other ranks
  void fetchFiberData();

  //! send vmValues data from fiberData_ back to the fibers where it belongs to and set in the respective field variable
  void updateFiberData();

  //! solve the 0D problem (Hodgkin-Huxley reaction term), starting from startTime
  void compute0D(double startTime, double timeStepWidth, int nTimeSteps);

  //! solve the 1D problem (diffusion), starting from startTime
  void compute1D(double startTime, double timeStepWidth, int nTimeSteps, double prefactor);

  //! compute the 0D-1D problem with Strang splitting
  void computeMonodomain();

  PythonConfig specificSettings_;    ///< config for this object

  NestedSolversType nestedSolvers_;   ///< the nested solvers object that would normally solve the problem

  /** data to be exchanged for computation of a single fiber
   *  The data stored herein is used for local computation.
   */
  struct FiberData
  {
    std::vector<double> elementLengths;   ///< lengths of the 1D elements
    std::vector<double> vmValues;         ///< values of Vm
    int valuesLength;                     ///< number of vmValues
    global_no_t valuesOffset;             ///< number of vmValues in previous entries in fiberData_
    int fiberNoGlobal;                    ///< fiberNo as given in settings (value of additionalArgument)
    int motorUnitNo;                      ///< motor unit no.
  };

  std::vector<FiberPointBuffers> fiberPointBuffers_;    ///< computation buffers for the 0D problem

  std::string fiberDistributionFilename_;  ///< filename of the fiberDistributionFile, which contains motor unit numbers for fiber numbers
  std::string firingTimesFilename_;        ///< filename of the firingTimesFile, which contains time points of stimulation for each motor unit

  std::vector<std::vector<bool>> firingEvents_;   ///< if a motor unit firingEvents_[timeStepNo][motorUnitNo]
  std::vector<int> motorUnitNo_;                  ///< number of motor unit for given fiber no motorUnitNo_[fiberNo]
  double setSpecificStatesCallFrequency_;         ///< value of option with the same name in the python settings
  double setSpecificStatesRepeatAfterFirstCall_;  ///< how long in ms the prescribed value should be set

  std::vector<FiberData> fiberData_;  ///< vector of fibers,
  int nFibersToCompute_;              ///< number of fibers where own rank is involved (>= n.fibers that are computed by own rank)
  int nInstancesToCompute_;           ///< number of instances of the Hodgkin-Huxley problem to compute on this rank
  double currentTime_;                ///< the current time used for the output writer
};