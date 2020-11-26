#pragma once

#include <Python.h>  // has to be the first included header
#include "time_stepping_scheme/00_time_stepping_scheme.h"

#include "easylogging++.h"

namespace TimeSteppingScheme
{

/** This class simply calls the run() method of the solver repeatedly.
 *  It is similar to a Coupling class, but only involves one sub-solver and no data transfer.
 */
template<typename Solver>
class RepeatedCallStatic :
  public TimeSteppingScheme
{
public:

  typedef typename Solver::SlotConnectorDataType SlotConnectorDataType;
  typedef typename Solver::Data Data;
  typedef typename Solver::FunctionSpace FunctionSpace;

  //! constructor
  RepeatedCallStatic(DihuContext context);

  //! advance simulation by the given time span [startTime_, endTime_] with given numberTimeSteps
  virtual void advanceTimeSpan();

  //! initialize solver
  void initialize();

  //! first initialize that run the stepping
  void run();

  //! return the data object, with the call to this method the output writers get the data to create their output files
  Data &data();

  //! Get the data that will be transferred in the operator splitting or coupling to the other term of the splitting/coupling.
  //! the transfer is done by the slot_connector_data_transfer class
  std::shared_ptr<SlotConnectorDataType> getSlotConnectorData();

private:

  Solver solver_;   //< the underlying solver object that will be stepped repeatedly
};

}  // namespace

#include "time_stepping_scheme/repeated_call_static.tpp"
