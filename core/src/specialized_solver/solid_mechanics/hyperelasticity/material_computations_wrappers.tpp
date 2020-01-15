#include "specialized_solver/solid_mechanics/hyperelasticity/hyperelasticity_solver.h"

#include <Python.h>  // has to be the first included header

namespace SpatialDiscretization
{

template<typename Term,int nDisplacementComponents>
void HyperelasticitySolver<Term,nDisplacementComponents>::
materialComputeInternalVirtualWork(
  std::shared_ptr<VecHyperelasticity> displacements,
  std::shared_ptr<VecHyperelasticity> internalVirtualWork
)
{
  // this is only for the static problem
  assert(nDisplacementComponents == 3);
  assert(displacements);
  assert(internalVirtualWork);

  // copy the values from the combined vector to displacements and pressure, this is needed for materialComputeResidual
  setDisplacementsAndPressureFromCombinedVec(displacements->valuesGlobal());

  // place the result vector to be used by materialComputeInternalVirtualWork
  Vec result = internalVirtualWork->valuesGlobal();
  PetscErrorCode ierr;
  ierr = VecSwap(result, solverVariableResidual_); CHKERRV(ierr);

  // compute the actual output of the nonlinear function
  materialComputeInternalVirtualWork();

  // reset the result vector
  VecSwap(result, solverVariableResidual_); CHKERRV(ierr);
}

template<typename Term,int nDisplacementComponents>
void HyperelasticitySolver<Term,nDisplacementComponents>::
solveForDisplacements(
  std::shared_ptr<VecHyperelasticity> externalVirtualWork,
  std::shared_ptr<VecHyperelasticity> displacements
)
{
  // solve ∂W_int - ∂W_ext = 0 with J = 1

  // assign new value for ∂W_ext, save previous vector
  Vec rhs = externalVirtualWork->valuesGlobal();
  PetscErrorCode ierr;
  ierr = VecSwap(rhs, externalVirtualWorkDead_); CHKERRV(ierr);

  Vec result = displacements->valuesGlobal();
  ierr = VecSwap(result, solverVariableSolution_); CHKERRV(ierr);

  nonlinearSolve();

  // restore value of ∂W_ext and solution
  ierr = VecSwap(rhs, externalVirtualWorkDead_); CHKERRV(ierr);
  ierr = VecSwap(result, solverVariableSolution_); CHKERRV(ierr);
}

template<typename Term,int nDisplacementComponents>
void HyperelasticitySolver<Term,nDisplacementComponents>::
solveDynamicProblem(
  std::shared_ptr<VecHyperelasticity> displacementsVelocitiesPressure, bool isFirstTimeStep,
  Vec internalVirtualWork, Vec &externalVirtualWorkDead, Vec accelerationTerm
)
{
  // solve δW_int - δW_ext,dead + ∫_Ω ρ0 1/dt (v^(n+1),L - v^(n),L)/dt ϕ^L ϕ^M δu dV = 0,
  // 1/dt (u^(n+1) - u^(n)) - v^(n+1) = 0
  // ∫ (J-1) δp = 0   ∀ δp,δu,δv

  assert(nDisplacementComponents == 6);

  PetscErrorCode ierr;
  Vec result = displacementsVelocitiesPressure->valuesGlobal();
  ierr = VecSwap(result, solverVariableSolution_); CHKERRV(ierr);

  // write reference output values
  if (isFirstTimeStep)
  {
    this->outputWriterManager_.writeOutput(this->data_, 0, 0.0);
    this->outputWriterManagerPressure_.writeOutput(this->pressureDataCopy_, 0, 0.0);
  }

  // save u,v,p of previous timestep
  setDisplacementsVelocitiesAndPressureFromCombinedVec(solverVariableSolution_,
                                                       this->data_.displacementsPreviousTimestep(),
                                                       this->data_.velocitiesPreviousTimestep(),
                                                       this->data_.pressurePreviousTimestep());

  LOG(DEBUG) << "saved initial uvp values: ";
  LOG(DEBUG) << *this->data_.displacementsPreviousTimestep();
  VLOG(1) << *this->data_.velocitiesPreviousTimestep();
  VLOG(1) << *this->data_.pressurePreviousTimestep();

  // find the solution for u,v,p of the nonlinear equation
  nonlinearSolve();

  // update geometry fields from displacements, compute PK2Stress, write output with output writers
  postprocessSolution();

  LOG(DEBUG) << "nonlinearSolve finished";

  // output with output writers of hyperelasticity_solver
  this->outputWriterManager_.writeOutput(this->data_, 0, endTime_);
  this->outputWriterManagerPressure_.writeOutput(this->pressureDataCopy_, 0, endTime_);

  // restore value of ∂W_ext and solution
  ierr = VecSwap(result, solverVariableSolution_); CHKERRV(ierr);

  // determine output values for δW
  externalVirtualWorkDead = externalVirtualWorkDead_;

  // compute internal virtual work
  setUVP(solverVariableSolution_);
  materialComputeInternalVirtualWork();
  ierr = VecCopy(solverVariableResidual_, internalVirtualWork); CHKERRV(ierr);

  // compute acceleration term
  ierr = VecZeroEntries(solverVariableResidual_); CHKERRV(ierr);
  materialAddAccelerationTermAndVelocityEquation();
  ierr = VecCopy(solverVariableResidual_, accelerationTerm); CHKERRV(ierr);

}

}  // namespace SpatialDiscretization
