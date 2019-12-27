#pragma once

#include <Python.h>  // has to be the first included header
#include <petscvec.h>

#include "interfaces/runnable.h"
#include "interfaces/splittable.h"
#include "interfaces/discretizable_in_time.h"
#include "control/dihu_context.h"
#include "control/python_config.h"

namespace Control
{

/** This class implements the time stepping, after each time step rebalance() is called, which can do load balancing.
 *  The actual implementation of the rebalancing has to be done in a derived class.
  */
template<class TimeStepping>
class LoadBalancingBase:
  public Runnable,
  public TimeSteppingScheme::TimeSteppingScheme
{
public:

  typedef typename TimeStepping::FunctionSpace FunctionSpace;
  typedef typename TimeStepping::Data Data;
  typedef typename TimeStepping::OutputConnectorDataType OutputConnectorDataType;

  //! constructor
  LoadBalancingBase(DihuContext context);

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
  OutputConnectorDataType &getOutputConnectorData();

protected:

  //! check if the degrees of freedom should be redistributed among the ranks
  virtual void rebalance() = 0;

  TimeStepping timeSteppingScheme_;   ///< the underlying timestepping method that is controlled by this class, e.g. Heun
};

}  // namespace

#include "control/load_balancing/load_balancing_base.tpp"
