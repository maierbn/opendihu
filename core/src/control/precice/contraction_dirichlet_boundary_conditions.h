#pragma once

#include <Python.h>  // has to be the first included header

#include "interfaces/runnable.h"
#include "data_management/specialized_solver/my_new_static_solver.h"   // adjust this include

#ifdef HAVE_PRECICE
#include "precice/SolverInterface.hpp"
#endif

namespace PreciceAdapter
{

/** Precice adapter for muscle.
 *
 *  The scheme is the following:
 *
 *  +--traction-(Neumann BC)-->[tendon]---displacement-(Dirichlet BC)-->[muscle]---stress--+
 *  |                                                                                      |
 *  +--------------------------------------------------------------------------------------+
 */
template<typename NestedSolver>
class ContractionDirichletBoundaryConditions :
  public Runnable
{
public:
  //! make the FunctionSpace of the NestedSolver class available
  typedef typename NestedSolver::TimeStepping2Type::FunctionSpace FunctionSpace;

  //! define the type of the data object,
  typedef typename NestedSolver::Data Data;   // either, if you do not need your own data object, use the data object of NestedSolver

  //! Define the type of data that will be transferred between solvers when there is a coupling scheme.
  //! Usually you define this type in the "Data" class and reuse it here.
  typedef typename NestedSolver::SlotConnectorDataType SlotConnectorDataType;

  //! constructor, gets the DihuContext object which contains all python settings
  ContractionDirichletBoundaryConditions(DihuContext context);

  //! initialize the object
  void initialize();

  //! run solution process, this first calls initialize() and then run()
  void run();

  //! reset state of this object, such that a new initialize() is necessary ("uninitialize")
  void reset();

  //! return the data object, with the call to this method the output writers get the data to create their output files
  Data &data();

  //! Get the data that will be transferred in the operator splitting or coupling to the other term of the splitting/coupling.
  //! the transfer is done by the slot_connector_data_transfer class
  std::shared_ptr<SlotConnectorDataType> getSlotConnectorData();

protected:

  //! initialize dirichlet boundary conditions of the hyperelasticity solver, set all Dirichlet boundary condition values that will be needed later to vector (0,0,0)
  void initializeDirichletBoundaryConditions();

#ifdef HAVE_PRECICE
  struct CouplingParticipant
  {
    bool isCouplingSurfaceBottom;                                      //< if the coupling surface where the tendon is attached is at the bottom of the muscle (z-)
    std::shared_ptr<precice::SolverInterface> preciceSolverInterface;  //< the precice solver interface that makes all preCICE functionality accessible
  };
  std::vector<CouplingParticipant> couplingParticipants_;     //< all coupling interface to different tendons

  //! read the data from the other partiticipant
  void preciceReadData(const CouplingParticipant &couplingParticipant);

  //! write the data to the other partiticipant
  void preciceWriteData(const CouplingParticipant &couplingParticipant);

#endif

  DihuContext context_;                       //< object that contains the python config for the current context and the global singletons meshManager and solverManager
  PythonConfig specificSettings_;             //< python object containing the value of the python config dict with corresponding key

  NestedSolver nestedSolver_;                 //< the nested solver that is controlled by this class
  std::shared_ptr<FunctionSpace> functionSpace_; //< 3D function space of displacements

  double maximumPreciceTimestepSize_;         //< maximum timestep size that precice will allow for the current time step
  double timeStepWidth_;                      //< timestep width of the solver
  int timeStepOutputInterval_;                //< interval in which to output current time

  bool haveCouplingSurfaceBottom_;            //< if there is a coupling surface at the bottom of the muscle
  bool haveCouplingSurfaceTop_;               //< if there is a coupling surface at the top of the muscle

  std::vector<int> preciceVertexIdsBottom_;   //< the vertex ids in precice of the geometry values at the bottom of the mesh
  std::vector<int> preciceVertexIdsTop_;      //< the vertex ids in precice of the geometry values at the top of the mesh
  int preciceMeshIdBottom_;                   //< mesh ID of precice of the surface mesh at the bottom that contains all nodes
  int preciceMeshIdTop_;                      //< mesh ID of precice of the surafce mesh at the top that contains all nodes
  int preciceDataIdDisplacements_;            //< data ID of precice of the displacements field to be exchanged
  int preciceDataIdVelocity_;                 //< data ID of precice of the velocity field to be exchanged
  int preciceDataIdTraction_;                 //< data ID of precice of the traction field to be exchanged

  int nNodesBottomSurfaceLocal_;              //< number of nodes on the bottom surface where the tendon is coupled
  int nNodesTopSurfaceLocal_;                 //< number of nodes on the top surface where the tendon is coupled

  std::vector<double> tractionValuesBottom_;  //< traction values to be transferred, in array-of-struct order (x,y,z,x,y,z,...)
  std::vector<double> tractionValuesTop_;     //< traction values to be transferred

  bool initialized_;                          //< if initialize() was already called
};

}  // namespace

#include "control/precice/contraction_dirichlet_boundary_conditions.tpp"
