// parallel-in-time (XBraid) Implicit Euler Solver ()

#pragma once

#include <Python.h>  // has to be the first included header

#include "interfaces/runnable.h"
#include "data_management/specialized_solver/PinT_MD.h"   // adjust this include
// #include "time_stepping_scheme/implicit_euler.h"
#include "specialized_solver/multidomain_solver/multidomain_solver.h"
#include "control/dihu_context.h"

#include "specialized_solver/parallel_in_time/MultiDomain/PinT_fun_MD.h"

namespace ParallelInTime {
template<class NestedSolverMD>
class PinTMD :
  public Runnable
{
public:
  //! make the FunctionSpace of the NestedSolverMD class available
  typedef typename NestedSolverMD::FunctionSpace FunctionSpace;

  //! define the type of the data object,
  //typedef typename NestedSolverMD::Data Data;   // either, if you do not need your own data object, use the data object of NestedSolverMD
  typedef ::Data::PinTMD<typename NestedSolverMD::FunctionSpace> Data;   // or, define your own data class, stored under "data_management/my_new_static_solver.h"

  //! Define the type of data that will be transferred between solvers when there is a coupling scheme.
  //! Usually you define this type in the "Data" class and reuse it here.
  typedef typename Data::OutputConnectorDataType OutputConnectorDataType;

  //! constructor, gets the DihuContext object which contains all python settings
  PinTMD(DihuContext context);

  //! initialize the object
  void initialize();

  //! run solution process, this first calls initialize() and then run()
  void run();

  //! reset state of this object, such that a new initialize() is necessary ("uninitialize")
  void reset();

  // void setTimeSpan(double tstart, double tstop);

  // void advanceTimeSpan();

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
      NestedSolverMD
    >
  > MultiDomainSolvers_;   //< vector of nested solvers (implicit euler) for solution on different grids

  std::vector<
    std::shared_ptr<Data>
  > data_;                                 //< the data objects that stores at least all field variables that should be output by output writers, one for each implicit euler

  std::shared_ptr<Partition::RankSubset> rankSubsetX_;   //< a rank subset containing the MPI communicator with numbering of ranks in space, created by xbraid

  bool initialized_;                          //< if initialize() was already called

  double tstart_ = 0.0;
  double tstop_ = 1.0;
  int ntime_ = 10;
  int nspace_ = 8;
  double xstart_ = 0.0;
  double xstop_ = 4.0;

  MPI_Comm communicatorTotal_= MPI_COMM_WORLD;
  MPI_Comm communicatorX_;
  MPI_Comm communicatorT_;

  // Braid variables
  braid_Core    core_;
  my_App       *app_;
  int print_level_   = 2;
  int max_levels_    = 4;
  int nrelax_ = 1;
  double tol_ = 1.0e-10;
  int fmg_ = 0;
  int cfactor_ = 2;
  int cfactor_first_ = -1;
  int nrelax_first_ = -1;

};
} // namespace ParallelInTime
#include "specialized_solver/parallel_in_time/MultiDomain/PinT_MD.tpp"
