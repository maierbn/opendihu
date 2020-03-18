#include "specialized_solver/fast_monodomain_solver/fast_monodomain_solver_base.h"

#include "partition/rank_subset.h"
#include "control/diagnostic_tool/stimulation_logging.h"

template<int nStates, int nIntermediates, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates,nIntermediates,DiffusionTimeSteppingScheme>::
reset()
{
  nestedSolvers_.reset();
}

template<int nStates, int nIntermediates, typename DiffusionTimeSteppingScheme>
typename FastMonodomainSolverBase<nStates,nIntermediates,DiffusionTimeSteppingScheme>::Data &FastMonodomainSolverBase<nStates,nIntermediates,DiffusionTimeSteppingScheme>::
data()
{
  return nestedSolvers_.data();
}

template<int nStates, int nIntermediates, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates,nIntermediates,DiffusionTimeSteppingScheme>::
setTimeSpan(double startTime, double endTime)
{
  nestedSolvers_.setTimeSpan(startTime, endTime);
}

//! get a reference to the nested solvers
template<int nStates, int nIntermediates, typename DiffusionTimeSteppingScheme>
typename FastMonodomainSolverBase<nStates,nIntermediates,DiffusionTimeSteppingScheme>::NestedSolversType &FastMonodomainSolverBase<nStates,nIntermediates,DiffusionTimeSteppingScheme>::
nestedSolvers()
{
  return nestedSolvers_;
}

template<int nStates, int nIntermediates, typename DiffusionTimeSteppingScheme>
std::shared_ptr<typename FastMonodomainSolverBase<nStates,nIntermediates,DiffusionTimeSteppingScheme>::OutputConnectorDataType>
FastMonodomainSolverBase<nStates,nIntermediates,DiffusionTimeSteppingScheme>::
getOutputConnectorData()
{
  return nestedSolvers_.getOutputConnectorData();
}
