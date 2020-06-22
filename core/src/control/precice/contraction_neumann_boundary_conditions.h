#pragma once

#include <Python.h>  // has to be the first included header

#include "interfaces/runnable.h"
#include "data_management/specialized_solver/my_new_static_solver.h"   // adjust this include

#ifdef HAVE_PRECICE
#include "precice/SolverInterface.hpp"
#endif

namespace PreciceAdapter
{

/** Precice adapter for tendon.
 *
 *  The scheme is the following:
 *
 *  +--traction-(Neumann BC)-->[tendon]---displacement-(Dirichlet BC)-->[muscle]---traction--+
 *  |                                                                                        |
 *  +----------------------------------------------------------------------------------------+
 */
template<typename NestedSolver>
class ContractionNeumannBoundaryConditions :
  public Runnable
{
public:
  //! make the FunctionSpace of the NestedSolver class available
  typedef typename NestedSolver::DisplacementsFunctionSpace FunctionSpace;

  //! define the type of the data object,
  typedef typename NestedSolver::Data Data;   // either, if you do not need your own data object, use the data object of NestedSolver

  //! constructor, gets the DihuContext object which contains all python settings
  ContractionNeumannBoundaryConditions(DihuContext context);

  //! initialize the object
  void initialize();

  //! run solution process, this first calls initialize() and then run()
  void run();

  //! reset state of this object, such that a new initialize() is necessary ("uninitialize")
  void reset();

  //! return the data object, with the call to this method the output writers get the data to create their output files
  Data &data();

protected:

#ifdef HAVE_PRECICE
  //! read the data from the other partiticipant
  void preciceReadData();

  //! write the data to the other partiticipant
  void preciceWriteData();

  std::unique_ptr<precice::SolverInterface> preciceSolverInterface_;  //< the precice solver interface that makes all preCICE functionality accessible
#endif

  DihuContext context_;                       //< object that contains the python config for the current context and the global singletons meshManager and solverManager
  PythonConfig specificSettings_;             //< python object containing the value of the python config dict with corresponding key

  NestedSolver nestedSolver_;                 //< the nested solver that is controlled by this class
  std::shared_ptr<FunctionSpace> functionSpace_;  //< the quadratic displacements function space of the dynamic solver

  double maximumPreciceTimestepSize_;         //< maximum timestep size that precice will allow for the current time step
  double timeStepWidth_;                      //< timestep width of the solver

  int outputConnectorSlotIdGamma_;            //< the number of the output connector slot that is used for gamma

  std::vector<int> preciceVertexIds_;         //< the vertex ids in precice of the geometry values
  int preciceMeshId_;                         //< mesh ID of precice of the mesh that contains all fiber nodes

  int preciceDataIdDisplacements_;            //< data ID of precice of the displacements field to be exchanged
  int preciceDataIdVelocity_;                 //< data ID of precice of the velocity field to be exchanged
  int preciceDataIdTraction_;                 //< data ID of precice of the traction field to be exchanged

  int nNodesSurfaceLocal_;                    //< number of nodes of the coupling surface
  bool isCouplingSurfaceBottom_;              //< the tendon is coupled at its bottom face to the muscle
  bool initialized_;                          //< if initialize() was already called
};

}  // namespace

#include "control/precice/contraction_neumann_boundary_conditions.tpp"
