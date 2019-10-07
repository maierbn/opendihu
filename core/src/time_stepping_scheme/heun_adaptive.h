#pragma once

#include "time_stepping_scheme/time_stepping_explicit.h"
#include "interfaces/runnable.h"
#include "data_management/time_stepping/time_stepping_heun.h"
#include "control/dihu_context.h"
#include <iostream>
#include <fstream>

namespace TimeSteppingScheme
{

/** The adaptive method to perform calculation with Heun method
 */
template<typename DiscretizableInTime>
class HeunAdaptive:
  public TimeSteppingExplicit<DiscretizableInTime>, public Runnable
{
public:

  //! constructor
  HeunAdaptive(DihuContext context);

  //! advance simulation by the given time span [startTime_, endTime_] with given numberTimeSteps, data in solution is used, afterwards new data is in solution
  void advanceTimeSpan();

  //! Returns the current time passed in the simulation. Used to trigger rebalancing
  double currentHeunTime();

private:
  // allowed tolerance
  double tolerance_;

  // offset to prevent small steps towards the end of the timeSpan
  double delta_;

  // safes timeStepWidth when it has to be cut because of the end of the timeSpan
  double savedTimeStepWidth_;

  // minimal timeStepWidth to use
  double minTimeStepWidth_;

  // check to show usage of minimal timeStepWidth
  bool minimum_check_;

  // check to show usage of delta
  bool ten_percent_check_;

  // value to adapt timeStepWidth with
  double alpha_;

  // estimator to estimate the error
  double estimator_;

  // lowest multiplier for modified method
  int lowestMultiplier_;

  // multiplicator for modified method
  int multiplicator_;

  // vector norm between two solutions
  PetscReal vecnorm;

  // temporary vectors for calculation
  Vec temp_solution_normal;
  Vec temp_solution_tilde;
  Vec temp_solution_tilde_intermediate;
  Vec temp_increment_1;
  Vec temp_increment_2;

  // counter for the successful steps
  int stepCounterSuccess_;

  // counter for the failed steps
  int stepCounterFail_;

  // read in performed adaption method
  std::string timeStepAdaptOption_;

  // current time of the simulation
  double currentTimeHeun_;

};
}  // namespace

#include "time_stepping_scheme/heun_adaptive.tpp"
