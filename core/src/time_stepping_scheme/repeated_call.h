#pragma once

#include <Python.h>  // has to be the first included header
#include "time_stepping_scheme/00_time_stepping_scheme.h"

#include "easylogging++.h"

namespace TimeSteppingScheme
{

/** This class simply calls the advanceTimeSpan() method of the solver repeatedly.
 *  It is similar to a Coupling class, but only involves one sub-solver and no data transfer.
 */
template<typename Solver>
class RepeatedCall :
  public TimeSteppingScheme
{
public:

  //! constructor
  RepeatedCall(DihuContext context);

  //! advance simulation by the given time span [startTime_, endTime_] with given numberTimeSteps
  virtual void advanceTimeSpan(bool withOutputWritersEnabled = true);

  //! initialize solver
  void initialize();

  //! first initialize that run the stepping
  void run();

  //! call the output writer on the data object, output files will contain currentTime, with callCountIncrement !=1 output timesteps can be skipped
  void callOutputWriter(int timeStepNo, double currentTime, int callCountIncrement = 1);

private:

  Solver solver_;   //< the underlying solver object that will be stepped repeatedly
};

}  // namespace

#include "time_stepping_scheme/repeated_call.tpp"
