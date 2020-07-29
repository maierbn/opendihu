#include "specialized_solver/fast_monodomain_solver/fast_monodomain_solver_base.h"

#include "partition/rank_subset.h"
#include "control/diagnostic_tool/stimulation_logging.h"

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
