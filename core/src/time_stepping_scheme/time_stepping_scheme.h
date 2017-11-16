#pragma once

#include "control/dihu_context.h"

namespace TimeSteppingScheme
{

class TimeSteppingScheme
{
public:
  TimeSteppingScheme(const DihuContext &context); 
 
  ///! advance simulation by the given time span [startTime_, endTime_] with given numberTimeSteps
  virtual void advanceTimeSpan() = 0;
 
  ///! set a new time step width, gets transferred to numberTimeSteps_
  void setTimeStepWidth(double timeStepWidth);
  
  ///! set a new number of time steps
  void setNumberTimeSteps(int numberTimeSteps);
  
  ///! set a new time interval that will be simulated by next call to advanceTimeSpan. This also potentially changes the time step width (it preserves the number of timesteps in the new time span)
  void setTimeSpan(double startTime, double endTime);

  ///! initialize time span from specificSettings_
  void initialize();
  
  ///! returns the Petsc solution vector
  virtual Vec &solution() = 0;
  
  ///! destructor  
  virtual ~TimeSteppingScheme() {}

protected:
 
  const DihuContext &context_;    ///< the context object containing everything to be stored
  
  double startTime_;        ///< start time of time interval to be simulated by call to advance
  double endTime_;          ///< end time of simulation
  int numberTimeSteps_;     ///< number of time steps in simulation time
  
  PyObject *specificSettings_;    ///< python object containing the value of the python config dict with corresponding key
};

}  // namespace