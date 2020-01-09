#pragma once

#include <Python.h>  // has to be the first included header
#include <petscvec.h>

#include "interfaces/runnable.h"
#include "interfaces/splittable.h"
#include "interfaces/discretizable_in_time.h"
#include "control/dihu_context.h"
#include "control/python_config.h"

/**
  */
template<class TimeStepping>
class MyNewTimesteppingSolver:
  public Runnable,
  public TimeSteppingScheme::TimeSteppingScheme
{
public:
  //! make the FunctionSpace of the TimeStepping class available
  typedef typename TimeStepping::FunctionSpace FunctionSpace;

  //! define the type of the data object,
  //typedef typename TimeStepping::Data Data;   // either, if you do not need your own data object, use the data object of TimeStepping
  typedef Data::MyNewTimesteppingSolver Data;   // or, define your own data class, stored under "data_management/my_new_timestepping_solver.h"

  //! Define the type of data that will be transferred between solvers when there is a coupling scheme.
  //! Usually you define this type in the "Data" class and reuse it here.
  typedef typename Data::OutputConnectorDataType OutputConnectorDataType;

  //! constructor, gets the DihuContext object which contains all python settings
  MyNewTimesteppingSolver(DihuContext context);

  //! advance simulation by the given time span [startTime_, endTime_] (set by setTimeSpan(), take a look at time_stepping_scheme/00_time_stepping_scheme.h)
  void advanceTimeSpan();

  //! initialize time span from specificSettings_
  void initialize();

  //! run solution process, this first calls initialize() and then run()
  void run();

  //! reset state of this object, such that a new initialize() is necessary
  void reset();

  //! return the data object of the timestepping scheme, with the call to this method the output writers get the data to create their output files
  Data &data();

  //! Get the data that will be transferred in the operator splitting or coupling to the other term of the splitting/coupling.
  //! the transfer is done by the output_connector_data_transfer class
  std::shared_ptr<OutputConnectorDataType> getOutputConnectorData();

protected:

  //! Here you can define private methods
  void executeMyHelperMethod();

  TimeStepping timeSteppingScheme_;   ///< the underlying timestepping method that is controlled by this class, e.g. Heun
};

#include "specialized_solver/my_new_solver/my_new_timestepping_solver.tpp"
