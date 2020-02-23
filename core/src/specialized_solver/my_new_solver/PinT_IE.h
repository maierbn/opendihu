#pragma once

#include <Python.h>  // has to be the first included header

#include "interfaces/runnable.h"
#include "data_management/specialized_solver/PinT_IE.h"   // adjust this include
#include "time_stepping_scheme/implicit_euler.h"
#include "control/dihu_context.h"

/** This is a template class that developers can copy and adjust to create their own solver.
 *  This solver is static, i.e. has no timestepping. For a dynamic solver, refer to the other template, "my_new_timestepping_solver.h".
 *  There are also the files "data_management/my_new_static_solver.{h,cpp}" that need to be adjusted.
 *  At the end, add an #include to this file in "opendihu.h".
 *
 *  Briefly explain what your solver does in this comment section.
  */
template<class NestedSolver>
class PinT :
  public Runnable
{
public:
  //! make the FunctionSpace of the NestedSolver class available
  typedef typename NestedSolver::FunctionSpace FunctionSpace;

  //! define the type of the data object,
  //typedef typename NestedSolver::Data Data;   // either, if you do not need your own data object, use the data object of NestedSolver
  typedef ::Data::PinT<typename NestedSolver::FunctionSpace> Data;   // or, define your own data class, stored under "data_management/my_new_static_solver.h"

  //! Define the type of data that will be transferred between solvers when there is a coupling scheme.
  //! Usually you define this type in the "Data" class and reuse it here.
  typedef typename Data::OutputConnectorDataType OutputConnectorDataType;

  //! constructor, gets the DihuContext object which contains all python settings
  PinT(DihuContext context);

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

  //! Here you can define private methods
  void PinT_initialize();
  void executeMyHelperMethod();

  DihuContext context_;                       //< object that contains the python config for the current context and the global singletons meshManager and solverManager
  OutputWriter::Manager outputWriterManager_; //< manager object holding all output writers
  PythonConfig specificSettings_;             //< python object containing the value of the python config dict with corresponding key


  std::vector<
    std::shared_ptr<
      NestedSolver
    >
  > implicitEulerSolvers_;   //< vector of nested solvers (implicit euler) for solution on different grids

  std::vector<
    std::shared_ptr<Data>
  > data_;                                 //< the data objects that stores at least all field variables that should be output by output writers, one for each implicit euler

  bool initialized_;                          //< if initialize() was already called

  double tstart_;
  double tstop_;
  int ntime_;
};

#include "specialized_solver/my_new_solver/PinT_IE.tpp"
