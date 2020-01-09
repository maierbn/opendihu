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
class MyNewStaticSolver:
  public Runnable,
  public TimeSteppingScheme::TimeSteppingScheme
{
public:

  typedef typename TimeStepping::FunctionSpace FunctionSpace;
  typedef typename TimeStepping::Data Data;
  typedef typename TimeStepping::OutputConnectorDataType OutputConnectorDataType;

  //! constructor
  MyNewStaticSolver(DihuContext context);

  //! advance simulation by the given time span [startTime_, endTime_]
  void advanceTimeSpan();

  //! initialize time span from specificSettings_
  void initialize();

  //! run solution process
  void run();

  //! reset state
  void reset();

  //! return the data object of the timestepping scheme
  Data &data();

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the output_connector_data_transfer class
  std::shared_ptr<OutputConnectorDataType> getOutputConnectorData();

protected:

  //! check if the degrees of freedom should be redistributed among the ranks
  virtual void rebalance() = 0;

  TimeStepping timeSteppingScheme_;   ///< the underlying timestepping method that is controlled by this class, e.g. Heun
};

#include "specialized_solver/my_new_solver/my_new_static_solver.tpp"
