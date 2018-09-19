#include "spatial_discretization/finite_element_method/05_time_stepping.h"
#include "discretizable_in_time/discretizable_in_time.h"

#include <Python.h>
#include <iostream>
#include <petscksp.h>
#include <vector>
#include <numeric>

#include "easylogging++.h"

#include "control/types.h"
#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "solver/solver_manager.h"
#include "solver/linear.h"

namespace SpatialDiscretization
{

template<typename FunctionSpaceType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType, Term>::
computeInverseMassMatrixTimesRightHandSide(Vec &result)
{
  // massMatrix * f_strong = rhs_weak
  Vec &rightHandSide = this->data_.rightHandSide()->valuesGlobal();   // rhs in weak formulation
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> massMatrix = this->data_.massMatrix();

  PetscErrorCode ierr;

  // create linear solver context
  VLOG(1) << "computeInverseMassMatrixTimesRightHandSide";
  initializeLinearSolver();

  // set matrix used for linear system and preconditioner to ksp context
  ierr = KSPSetOperators(*ksp_, massMatrix->valuesGlobal(), massMatrix->valuesGlobal()); CHKERRV(ierr);

#if 0
  // non zero initial values
  PetscScalar scalar = .5;
  ierr = VecSet(result, scalar); CHKERRV(ierr);
  ierr = KSPSetInitialGuessNonzero(*ksp_, PETSC_TRUE); CHKERRV(ierr);
#endif

  // solve the system, KSP assumes the initial guess is to be zero (and thus zeros it out before solving)
  ierr = KSPSolve(*ksp_, rightHandSide, result); CHKERRV(ierr);

  int numberOfIterations = 0;
  PetscReal residualNorm = 0.0;
  ierr = KSPGetIterationNumber(*ksp_, &numberOfIterations); CHKERRV(ierr);
  ierr = KSPGetResidualNorm(*ksp_, &residualNorm); CHKERRV(ierr);

  KSPConvergedReason convergedReason;
  ierr = KSPGetConvergedReason(*ksp_, &convergedReason); CHKERRV(ierr);

  VLOG(1) << "Rhs (" << this->data_.rightHandSide()->nDofsGlobal() << " global dofs) recovered in " << numberOfIterations << " iterations, residual norm " << residualNorm
    << ": " << PetscUtility::getStringLinearConvergedReason(convergedReason);
}

template<typename FunctionSpaceType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<FunctionSpaceType, QuadratureType, Term>::
evaluateTimesteppingRightHandSideExplicit(Vec &input, Vec &output, int timeStepNo, double currentTime)
{
  // this method computes output = M^{-1}*K*input
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> stiffnessMatrix = this->data_.stiffnessMatrix();
  Vec &rhs = this->data_.rightHandSide()->valuesGlobal();

  // check if matrix and vector sizes match
  PetscUtility::checkDimensionsMatrixVector(stiffnessMatrix->valuesGlobal(), input);

  // compute rhs = stiffnessMatrix*input
  PetscErrorCode ierr;
  ierr = MatMult(stiffnessMatrix->valuesGlobal(), input, rhs); CHKERRV(ierr);

  // compute output = massMatrix^{-1}*rhs
  computeInverseMassMatrixTimesRightHandSide(output);

  // output data as specified in outputWriter
  this->data_.print();
  this->outputWriterManager_.writeOutput(this->data_, timeStepNo, currentTime);
}

} // namespace SpatialDiscretization
