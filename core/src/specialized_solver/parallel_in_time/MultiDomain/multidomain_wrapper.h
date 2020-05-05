#pragma once

#include <Python.h>  // has to be the first included header

#include "interfaces/runnable.h"
#include "time_stepping_scheme/00_time_stepping_scheme.h"
#include "data_management/specialized_solver/PinT_MD.h"   // adjust this include

/** A wrapper around the multidomain_no_fat scheme to be used within the parallel in time class.
 *  A valid class for StringSplittingMultidomain would be:
 * ```
 *   OperatorSplitting::Strang<
 *     Control::MultipleInstances<
 *       TimeSteppingScheme::Heun<
 *         CellmlAdapter<
 *           4,9,  // nStates,nIntermediates: 57,1 = Shorten, 4,9 = Hodgkin Huxley
 *           FunctionSpace::FunctionSpace<MeshType,BasisFunction::LagrangeOfOrder<1>>  // same function space as for anisotropic diffusion
 *         >
 *       >
 *     >,
 *     TimeSteppingScheme::MultidomainSolver<              // multidomain
 *       SpatialDiscretization::FiniteElementMethod<       //FEM for initial potential flow, fibre directions
 *         MeshType,
 *         BasisFunction::LagrangeOfOrder<1>,
 *         Quadrature::Gauss<3>,
 *         Equation::Static::Laplace
 *       >,
 *       SpatialDiscretization::FiniteElementMethod<   // anisotropic diffusion
 *         MeshType,
 *         BasisFunction::LagrangeOfOrder<1>,
 *         Quadrature::Gauss<5>,
 *         Equation::Dynamic::DirectionalDiffusion
 *       >
 *     >
 *   >
 * ```
  */
template<class StrangSplittingMultidomain>
class MultidomainWrapper :
  public Runnable
{
public:
  //! make the FunctionSpace of the StrangSplittingMultidomain class available
  typedef typename StrangSplittingMultidomain::FunctionSpace FunctionSpace;

  //! define the type of the data object,
  typedef typename StrangSplittingMultidomain::Data Data;   // data object

  //! Define the type of data that will be transferred between solvers when there is a coupling scheme.
  //! Usually you define this type in the "Data" class and reuse it here.
  typedef typename StrangSplittingMultidomain::OutputConnectorDataType OutputConnectorDataType;

  //! constructor, gets the DihuContext object which contains all python settings
  MultidomainWrapper(DihuContext context);

  //! advance simulation by the given time span [startTime_, endTime_] (set by setTimeSpan(), take a look at time_stepping_scheme/00_time_stepping_scheme.h)
  void advanceTimeSpan();

  //! initialize time span from specificSettings_
  void initialize();

  //! run solution process, this first calls initialize() and then run()
  void run();

  //! setup a new system matrix within the multidomain solver
  void setSystemMatrix(double timeStepWidth);

  //! get the local number of solution values, which is the number of entries set by getSolution
  int nSolutionValuesLocal();

  //! get a copy of the current solution as a Petsc Vec that contains all entries, i.e. all states for all multidomain compartments
  //! The layout is such that state no. increases fastest, then compartment no., then value no.
  //! [compartment0state0value0, compartment0state1value0, ..., compartment0stateNvalue0, ..., compartemen1state0value0, ..., compartemen1stateNvalue0, ..., compartment0state0value1, compartement0state1value1, ...]
  void getSolution(double *data);

  //! copy the values given in solution back to the internal solution variables (states)
  //! @param solution This Vec should be the one that was returned by getSolutionAsVec.
  void setSolution(double *data);

  //! set a new time step width, gets transferred to numberTimeSteps_
  void setTimeStepWidth(double timeStepWidth);

  //! set a new number of time steps
  void setNumberTimeSteps(int numberTimeSteps);

  //! set a new time interval that will be simulated by next call to advanceTimeSpan. This also potentially changes the time step width (it preserves the number of timesteps in the new time span)
  void setTimeSpan(double startTime, double endTime);

  //! interval for output of time step number and time
  int timeStepOutputInterval();

  //! start time of time interval to be simulated
  double startTime();

  //! end time of simulation
  double endTime();

  //! number of time steps in simulation time
  int numberTimeSteps();

  //! time step for simulation
  double timeStepWidth();

  //! reset state of this object, such that a new initialize() is necessary ("uninitialize")
  void reset();

  //! return the data object of the timestepping scheme, with the call to this method the output writers get the data to create their output files
  Data &data();

  //! Get the data that will be transferred in the operator splitting or coupling to the other term of the splitting/coupling.
  //! the transfer is done by the output_connector_data_transfer class
  std::shared_ptr<OutputConnectorDataType> getOutputConnectorData();

protected:

  StrangSplittingMultidomain strangSplittingMultidomain_;   //< the underlying strang splitting scheme with multidomain
};

#include "specialized_solver/parallel_in_time/MultiDomain/multidomain_wrapper.tpp"
