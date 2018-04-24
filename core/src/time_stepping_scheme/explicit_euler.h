#pragma once

#include <Python.h>  // has to be the first included header
#include "time_stepping_scheme/time_stepping_scheme_ode.h"
#include "control/runnable.h"
#include "data_management/time_stepping.h"
#include "control/dihu_context.h"
#include "model_order_reduction/pod.h"

namespace TimeSteppingScheme
{

/** The explicit Euler integration scheme (forward Euler), u_{t+1} = u_{t} + dt*f(t)
  */
template<typename DiscretizableInTime>
class ExplicitEuler : 
  public TimeSteppingSchemeOde<DiscretizableInTime>, public Runnable
{
public:
 
  //! constructor
  ExplicitEuler(DihuContext context);
  
  //! advance simulation by the given time span [startTime_, endTime_] with given numberTimeSteps, data in solution is used, afterwards new data is in solution
  void advanceTimeSpan();
  
  //! run the simulation
  void run();
private: 
};

}  // namespace

#include "time_stepping_scheme/explicit_euler.tpp"
