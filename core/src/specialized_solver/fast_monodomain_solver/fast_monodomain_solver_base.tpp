#include "specialized_solver/fast_monodomain_solver/fast_monodomain_solver_base.h"

#include "partition/rank_subset.h"
#include "control/diagnostic_tool/stimulation_logging.h"

template<int nStates, int nIntermediates>
void FastMonodomainSolverBase<nStates,nIntermediates>::
reset()
{
  nestedSolvers_.reset();
}

template<int nStates, int nIntermediates>
typename FastMonodomainSolverBase<nStates,nIntermediates>::Data &FastMonodomainSolverBase<nStates,nIntermediates>::
data()
{
  return nestedSolvers_.data();
}

template<int nStates, int nIntermediates>
void FastMonodomainSolverBase<nStates,nIntermediates>::
setTimeSpan(double startTime, double endTime)
{
  nestedSolvers_.setTimeSpan(startTime, endTime);
}

//! get a reference to the nested solvers
template<int nStates, int nIntermediates>
typename FastMonodomainSolverBase<nStates,nIntermediates>::NestedSolversType &FastMonodomainSolverBase<nStates,nIntermediates>::
nestedSolvers()
{
  return nestedSolvers_;
}

template<int nStates, int nIntermediates>
std::shared_ptr<typename FastMonodomainSolverBase<nStates,nIntermediates>::OutputConnectorDataType>
FastMonodomainSolverBase<nStates,nIntermediates>::
getOutputConnectorData()
{
  return nestedSolvers_.getOutputConnectorData();
}
