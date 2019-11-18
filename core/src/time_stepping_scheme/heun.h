#pragma once

#include "time_stepping_scheme/time_stepping_explicit.h"
#include "interfaces/runnable.h"
#include "data_management/time_stepping/time_stepping_heun.h"
#include "control/dihu_context.h"

namespace TimeSteppingScheme
{

/** The Heun integration scheme, u_{t+1} = u_{t} + 0.5*dt*(f(u_{t})+f(u*))
 *
 *  However, we compute it in the way: u_{t+1} = u* + 0.5*dt*(f(u*)-f(u_{t})) where u* = u_{t} + dt*f(u_{t}).
 *  (more round off this way, but less storage required)
 *
 */
template<typename DiscretizableInTime>
class Heun :
  public TimeSteppingExplicit<DiscretizableInTime>, public Runnable
{
public:

  //! constructor
  Heun(DihuContext context);

  //! initialize the data object
  virtual void initialize();

  //! advance simulation by the given time span [startTime_, endTime_] with given numberTimeSteps, data in solution is used, afterwards new data is in solution
  void advanceTimeSpan();

  //! run the simulation
  void run();

};

}  // namespace

#include "time_stepping_scheme/heun.tpp"
