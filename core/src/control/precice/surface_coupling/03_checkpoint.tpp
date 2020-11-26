#include "control/precice/surface_coupling/03_checkpoint.h"

#include <sstream>

namespace Control
{

//! save the current state of the simulation to a checkpoint
template<typename NestedSolver>
void PreciceAdapterCheckpoint<NestedSolver>::
saveCheckpoint(double currentTime)
{
  savedCurrentTime_ = currentTime;
  LOG(DEBUG) << "save checkpoint at time " << currentTime;

  PetscErrorCode ierr;

  // if the variable savedState_ does not yet exist, create as duplicate from the state vector
  if (savedState_ == PETSC_NULL)
  {
    ierr = VecDuplicate(this->currentState(this->nestedSolver_), &savedState_); CHKERRV(ierr);
  }

  // copy values
  ierr = VecCopy(this->currentState(this->nestedSolver_), savedState_); CHKERRV(ierr);

}

//! load all values from the last saved checkpoint
//! @return current simulation time
template<typename NestedSolver>
double PreciceAdapterCheckpoint<NestedSolver>::
loadCheckpoint()
{
  LOG(DEBUG) << "load checkpoint from time " << savedCurrentTime_;

  PetscErrorCode ierr;

  // copy values back
  ierr = VecCopy(savedState_, this->currentState(this->nestedSolver_));
  if (ierr != 0)
    LOG(FATAL) << "Loading checkpoint failed.";


  return savedCurrentTime_;
}

}  // namespace
