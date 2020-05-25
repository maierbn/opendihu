#pragma once

#include <Python.h>  // has to be the first included header

#include "interfaces/runnable.h"
#include "data_management/specialized_solver/my_new_static_solver.h"   // adjust this include

#ifdef HAVE_PRECICE
#include "precice/SolverInterface.hpp"
#endif

namespace PreciceAdapter
{

/** Precice adapter for partitioned fibers.
  */
template<class NestedSolver>
class ContractionNeumannBoundaryConditions :
  public Runnable
{
public:
  //! make the FunctionSpace of the NestedSolver class available
  typedef typename NestedSolver::FunctionSpace FunctionSpace;

  //! define the type of the data object,
  typedef typename NestedSolver::Data Data;   // either, if you do not need your own data object, use the data object of NestedSolver

  //! Define the type of data that will be transferred between solvers when there is a coupling scheme.
  //! Usually you define this type in the "Data" class and reuse it here.
  typedef typename NestedSolver::OutputConnectorDataType OutputConnectorDataType;

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

  //! Get the data that will be transferred in the operator splitting or coupling to the other term of the splitting/coupling.
  //! the transfer is done by the output_connector_data_transfer class
  std::shared_ptr<OutputConnectorDataType> getOutputConnectorData();

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

  double maximumPreciceTimestepSize_;         //< maximum timestep size that precice will allow for the current time step
  double timeStepWidth_;                      //< timestep width of the solver

  int outputConnectorSlotIdGamma_;            //< the number of the output connector slot that is used for gamma

  std::vector<int> preciceVertexIds_;         //< the vertex ids in precice of the geometry values
  int preciceMeshId_;                         //< mesh ID of precice of the mesh that contains all fiber nodes

  int preciceDataIdGeometry_;                 //< data ID of precice of the geometry information to be exchanged
  int preciceDataIdGamma_;                    //< data ID of precice of the gamma field to be exchanged
  int preciceDataIdLambda_;                   //< data ID of precice of the lambda field to be exchanged
  int preciceDataIdLambdaDot_;                //< data ID of precice of the lambda dot field to be exchanged

  bool initialized_;                          //< if initialize() was already called
};

}  // namespace

#include "control/precice/muscle_contraction.tpp"
