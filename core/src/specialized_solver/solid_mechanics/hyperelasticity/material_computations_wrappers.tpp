#include "specialized_solver/solid_mechanics/hyperelasticity/hyperelasticity_solver.h"

#include <Python.h>  // has to be the first included header

namespace SpatialDiscretization
{

template<typename Term,int nDisplacementComponents>
void HyperelasticitySolver<Term,nDisplacementComponents>::
materialComputeInternalVirtualWork(
  std::shared_ptr<PetscVec> displacements,
  std::shared_ptr<PetscVec> internalVirtualWork
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
  std::shared_ptr<PetscVec> externalVirtualWork,
  std::shared_ptr<PetscVec> displacements
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
  std::shared_ptr<PetscVec> displacementsVelocitiesPressure
)
{
  // solve δW_int - δW_ext,dead + ∫_Ω ρ0 1/dt (v^(n+1),L - v^(n),L)/dt ϕ^L ϕ^M δu dV = 0,
  // 1/dt (u^(n+1) - u^(n)) - v^(n+1) = 0
  // ∫ (J-1) δp = 0   ∀ δp,δu,δv

  assert(nDisplacementComponents == 6);

  PetscErrorCode ierr;
  Vec result = displacementsVelocitiesPressure->valuesGlobal();
  ierr = VecSwap(result, solverVariableSolution_); CHKERRV(ierr);

  // save u,v,p of previous timestep
  setDisplacementsVelocitiesAndPressureFromCombinedVec(result,
                                                       this->data_.displacementsPreviousTimestep(),
                                                       this->data_.velocitiesPreviousTimestep()
                                                       this->data_.pressurePreviousTimestep());

  // find the solution for u,v,p of the nonlinear equation
  nonlinearSolve();

  // restore value of ∂W_ext and solution
  ierr = VecSwap(result, solverVariableSolution_); CHKERRV(ierr);
}

}  // namespace SpatialDiscretization
