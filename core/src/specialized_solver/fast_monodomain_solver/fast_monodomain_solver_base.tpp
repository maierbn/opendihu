#include "specialized_solver/fast_monodomain_solver/fast_monodomain_solver_base.h"

#include "partition/rank_subset.h"
#include "control/diagnostic_tool/stimulation_logging.h"

//! call the output writer on the data object, output files will contain currentTime, with callCountIncrement !=1 output timesteps can be skipped
template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
callOutputWriter(int timeStepNo, double currentTime, int callCountIncrement)
{
  // call output writer of diffusion
  std::vector<typename NestedSolversType::TimeSteppingSchemeType> &instances = nestedSolvers_.instancesLocal();

  for (int i = 0; i < instances.size(); i++)
  {
    // call write output of MultipleInstances, callCountIncrement is the number of times the output writer would have been called without FastMonodomainSolver
    instances[i].timeStepping2().writeOwnOutput(0, currentTime_, nTimeStepsSplitting_);
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

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates,nAlgebraics,DiffusionTimeSteppingScheme>::
loadFiberState()
{


  //fetchFiberData();

  // LOG(INFO) << "current vm values";
  // for (int i=0; i<fiberData_.size(); i++){
  //   LOG(INFO) << "fiber " << i << " has vm " << fiberData_[i].vmValues;
  // }

  // LOG(INFO) << "saved vm values";
  // for (int i=0; i<fiberDataOld_.size(); i++){
  //   LOG(INFO) << "fiber " << i << " has vm " << fiberDataOld_[i].vmValues;
  // }

  // LOG(INFO) << "update fiber state";
  fiberData_ = fiberDataOld_;
  fiberHasBeenStimulated_ = fiberHasBeenStimulatedOld_;
  firingEvents_ = firingEventsOld_;
  // for (int i=0; i<fiberData_.size(); i++){
  //   LOG(INFO) << "fiber " << i << " has new vm " << fiberData_[i].vmValues;
  // }
  fiberPointBuffersParameters_ =fiberPointBuffersParametersOld_ ;        //< constant parameter values, changing parameters is not implemented
  fiberPointBuffersAlgebraicsForTransfer_ = fiberPointBuffersAlgebraicsForTransferOld_;
  statesForTransferIndices_ = statesForTransferIndicesOld_;
  algebraicsForTransferIndices_ = algebraicsForTransferIndicesOld_;
  fiberPointBuffersStatesAreCloseToEquilibrium_= fiberPointBuffersStatesAreCloseToEquilibriumOld_;
  
  gpuParameters_ = gpuParametersOld_;              //< for "gpu": constant parameter values, in struct of array memory layout: gpuParameters_[parameterNo*nInstances + instanceNo]
  gpuAlgebraicsForTransfer_ = gpuAlgebraicsForTransferOld_;   //< for "gpu": algebraic values to use for slot connector data, in struct of array memory layout: gpuAlgebraicsForTransfer_[algebraicNo*nInstances + instanceNo]
  gpuStatesForTransfer_ = gpuStatesForTransferOld_;       //< for "gpu": state values to use for slot connector data, in struct of array memory layout: gpuStatesForTransfer_[stateInThisListIndex*nInstances + instanceNo]
  gpuElementLengths_ = gpuElementLengthsOld_;          //< for "gpu": the lengths of the 1D elements, in struct of array memory layout: gpuElementLengths_[fiberDataNo*nElementsOnFiber + elementNo]
  gpuFiringEvents_ = gpuFiringEventsOld_;              //< for "gpu": if a motor unit fires at a specified time, 1=yes, 0=no, gpuFiringEvents_[timeStepNo*nMotorUnits + motorUnitNo]
  gpuSetSpecificStatesFrequencyJitter_ = gpuSetSpecificStatesFrequencyJitterOld_;

  gpuFiberIsCurrentlyStimulated_ = gpuFiberIsCurrentlyStimulatedOld_; //< for "gpu": the value of fiberData_[].currentlyStimulating
  gpuMotorUnitNo_ =gpuMotorUnitNoOld_;                      //< motor unit no.
  gpuFiberStimulationPointIndex_ = gpuFiberStimulationPointIndexOld_;       //< index of the point on the fiber where to stimulate, i.e. position of the neuromuscular junction, if at center, it is equal to (int)(fiberData_[fiberDataNo].valuesLength / 2)
  gpuLastStimulationCheckTime_ = gpuLastStimulationCheckTimeOld_;      //< last time the fiber was checked for stimulation
  gpuSetSpecificStatesCallFrequency_ = gpuSetSpecificStatesCallFrequencyOld_;        //< value of option with the same name in the python settings
  gpuSetSpecificStatesRepeatAfterFirstCall_ = gpuSetSpecificStatesRepeatAfterFirstCallOld_; //< how long in ms the prescribed value should be set
  gpuSetSpecificStatesCallEnableBegin_ = gpuSetSpecificStatesCallEnableBeginOld_;      //< value of option with the same name in the python settings
  gpuCurrentJitter_ = gpuCurrentJitterOld_;                         //< current absolute value of jitter to add to setSpecificStatesCallFrequency
  gpuJitterIndex_ = gpuJitterIndexOld_;                              //< index of the vector in setSpecificStatesFrequencyJitter which is the current value to use
  gpuVmValues_ = gpuVmValuesOld_;     
   
  
  updateFiberData();
}

template<int nStates, int nAlgebraics, typename DiffusionTimeSteppingScheme>
void FastMonodomainSolverBase<nStates, nAlgebraics, DiffusionTimeSteppingScheme>::
saveFiberState()
{

  fetchFiberData();
  
  // LOG(INFO) << "current fiberHasBeenStimulated_ values" << fiberHasBeenStimulated_;
  // LOG(INFO) << "current firing events values" << firingEvents_;

  // LOG(INFO) << "save fiber state";
  fiberDataOld_ = fiberData_;
  fiberHasBeenStimulatedOld_ = fiberHasBeenStimulated_;
  //fiberHasBeenStimulatedOld_ = fiberHasBeenStimulated_;
  firingEventsOld_ = firingEvents_;
  fiberPointBuffersParametersOld_ =fiberPointBuffersParameters_ ;        //< constant parameter values, changing parameters is not implemented
  fiberPointBuffersAlgebraicsForTransferOld_ = fiberPointBuffersAlgebraicsForTransfer_;
  statesForTransferIndicesOld_ = statesForTransferIndices_;
  algebraicsForTransferIndicesOld_ = algebraicsForTransferIndices_;
  fiberPointBuffersStatesAreCloseToEquilibriumOld_ = fiberPointBuffersStatesAreCloseToEquilibrium_;  
  
  gpuParametersOld_ = gpuParameters_;              //< for "gpu": constant parameter values, in struct of array memory layout: gpuParameters_[parameterNo*nInstances + instanceNo]
  gpuAlgebraicsForTransferOld_ = gpuAlgebraicsForTransfer_;   //< for "gpu": algebraic values to use for slot connector data, in struct of array memory layout: gpuAlgebraicsForTransfer_[algebraicNo*nInstances + instanceNo]
  gpuStatesForTransferOld_ = gpuStatesForTransfer_;       //< for "gpu": state values to use for slot connector data, in struct of array memory layout: gpuStatesForTransfer_[stateInThisListIndex*nInstances + instanceNo]
  gpuElementLengthsOld_ = gpuElementLengths_;          //< for "gpu": the lengths of the 1D elements, in struct of array memory layout: gpuElementLengths_[fiberDataNo*nElementsOnFiber + elementNo]
  gpuFiringEventsOld_ = gpuFiringEvents_;              //< for "gpu": if a motor unit fires at a specified time, 1=yes, 0=no, gpuFiringEvents_[timeStepNo*nMotorUnits + motorUnitNo]
  gpuSetSpecificStatesFrequencyJitterOld_ = gpuSetSpecificStatesFrequencyJitter_;

  gpuFiberIsCurrentlyStimulatedOld_ = gpuFiberIsCurrentlyStimulated_; //< for "gpu": the value of fiberData_[].currentlyStimulating
  gpuMotorUnitNoOld_ = gpuMotorUnitNo_;                      //< motor unit no.
  gpuFiberStimulationPointIndexOld_ = gpuFiberStimulationPointIndex_;       //< index of the point on the fiber where to stimulate, i.e. position of the neuromuscular junction, if at center, it is equal to (int)(fiberData_[fiberDataNo].valuesLength / 2)
  gpuLastStimulationCheckTimeOld_ = gpuLastStimulationCheckTime_;      //< last time the fiber was checked for stimulation
  gpuSetSpecificStatesCallFrequencyOld_ = gpuSetSpecificStatesCallFrequency_;        //< value of option with the same name in the python settings
  gpuSetSpecificStatesRepeatAfterFirstCallOld_ = gpuSetSpecificStatesRepeatAfterFirstCall_; //< how long in ms the prescribed value should be set
  gpuSetSpecificStatesCallEnableBeginOld_ = gpuSetSpecificStatesCallEnableBegin_;      //< value of option with the same name in the python settings
  gpuCurrentJitterOld_ = gpuCurrentJitter_;                         //< current absolute value of jitter to add to setSpecificStatesCallFrequency
  gpuJitterIndexOld_ = gpuJitterIndex_;                              //< index of the vector in setSpecificStatesFrequencyJitter which is the current value to use
  gpuVmValuesOld_ = gpuVmValues_;     
   
   //<  [fiberPointNo][algebraicToTransferNo], algebraic values to use for slot connector data
  // LOG(INFO) << "saved vm values";
  // LOG(INFO) << "saved fiberHasBeenStimulated_ values" << fiberHasBeenStimulatedOld_;
  // LOG(INFO) << "saved firing events values" << firingEventsOld_;
  
  //updateFiberData();

}