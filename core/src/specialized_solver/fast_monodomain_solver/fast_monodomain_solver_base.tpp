#include "specialized_solver/fast_monodomain_solver/fast_monodomain_solver_base.h"

#include "partition/rank_subset.h"
#include "control/diagnostic_tool/stimulation_logging.h"

//! call the output writer on the data object, output files will contain currentTime, with callCountIncrement !=1 output timesteps can be skipped
template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
callOutputWriter(int timeStepNo, double currentTime, int callCountIncrement)
{
  // call output writer of diffusion
  LOG(DEBUG)<<"callOutputWriter with timeStepNo, currentTime = " << timeStepNo << " " << currentTime;

  std::vector<typename NestedSolversType::TimeSteppingSchemeType> &instances = nestedSolvers_.instancesLocal();

  for (int i = 0; i < instances.size(); i++)
  {
    // call write output of MultipleInstances, callCountIncrement is the number of times the output writer would have been called without FastMonodomainSolver
    instances[i].timeStepping2().writeOwnOutput(timeStepNo, currentTime_, nTimeStepsSplitting_);
  }
  
}

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
reset()
{
  nestedSolvers_.reset();
}

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
typename FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::Data &FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
data()
{
  return nestedSolvers_.data();
}

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
setTimeSpan(double startTime, double endTime)
{
  nestedSolvers_.setTimeSpan(startTime, endTime);
}

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
double FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
startTime() {
  return nestedSolvers_.startTime();
}

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
double FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
endTime() {
  return nestedSolvers_.endTime();
}

//! get a reference to the nested solvers
template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
typename FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::NestedSolversType &FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
nestedSolvers()
{
  return nestedSolvers_;
}

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
std::shared_ptr<typename FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::SlotConnectorDataType>
FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
getSlotConnectorData()
{
  return nestedSolvers_.getSlotConnectorData();
}

