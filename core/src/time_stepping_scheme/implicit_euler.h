#pragma once

#include <Python.h>  // has to be the first included header
#include "time_stepping_scheme/time_stepping_scheme_ode.h"
#include "control/runnable.h"
#include "data_management/time_stepping.h"
#include "control/dihu_context.h"

namespace TimeSteppingScheme
{

/** The implicit Euler integration scheme (backward Euler), u_{t+1} = u_{t} + dt*f(t+1)
  */
template<typename DiscretizableInTime>
class ImplicitEuler :
  public TimeSteppingSchemeOde<DiscretizableInTime>, public Runnable
{
public:

  //! constructor
  ImplicitEuler(DihuContext context);

  //! advance simulation by the given time span [startTime_, endTime_] with given numberTimeSteps, data in solution is used, afterwards new data is in solution
  void advanceTimeSpan();
  
  //! initialize discretizableInTime
  void initialize();

  //! run the simulation
  void run();
private:
  bool initialized_=false;     ///< if initialize() was already called
};

}  // namespace

#include "time_stepping_scheme/implicit_euler.tpp"
