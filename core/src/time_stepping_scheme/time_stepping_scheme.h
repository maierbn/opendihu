#pragma once

#include <Python.h>  // has to be the first included header
#include "control/dihu_context.h"
#include "output_writer/manager.h"
#include "interfaces/splitable.h"

#include "easylogging++.h"

namespace TimeSteppingScheme
{

class TimeSteppingScheme :
  public Splitable
{
public:

  //! constructor
  TimeSteppingScheme(DihuContext context);

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
  
  //! interval for output of time step number and time
  int timeStepOutputInterval();

  ///! reset state such that new initialization becomes necessary
  virtual void reset();

  ///! return whether the scheme has a specified mesh type and is not independent of the mesh type
  virtual bool knowsMeshType() = 0;

  ///! start time of time interval to be simulated
  double startTime();
  
  ///! end time of simulation
  double endTime();
  
  ///! number of time steps in simulation time
  double numberTimeSteps();
  
  ///! time step for simulation
  double timeStepWidth();
  
  ///! python object containing the value of the python config dict with corresponding key
  PyObject *specificSettings();

  ///! destructor
  virtual ~TimeSteppingScheme() {}

protected:

  DihuContext context_;    ///< object that contains the python config for the current context and the global singletons meshManager and solverManager
  OutputWriter::Manager outputWriterManager_; ///< manager object holding all output writer
  int timeStepOutputInterval_;    ///< time step number and time is output every timeStepOutputInterval_ time steps

  double startTime_;        ///< start time of time interval to be simulated by call to advance
  double endTime_;          ///< end time of simulation

  int numberTimeSteps_;     ///< number of time steps in simulation time, this one is always used to compute the real timeStepWidth

  bool isTimeStepWidthSignificant_;   ///< if time step width will be used to determine number of steps
  double timeStepWidth_;        ///< a timeStepWidth value that is used to compute the number of time steps
  std::string durationLogKey_;   ///< the key under which the duration of the time stepping is saved in the log

  PyObject *specificSettings_;    ///< python object containing the value of the python config dict with corresponding key
  bool initialized_;      ///< if initialize() was already called
};

}  // namespace
