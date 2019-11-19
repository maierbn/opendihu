#include "specialized_solver/solid_mechanics/hyperelasticity/hyperelasticity_solver.h"

#include <Python.h>  // has to be the first included header

namespace SpatialDiscretization
{

template<typename Term>
void HyperelasticitySolver<Term>::
materialComputeInternalVirtualWork(
  std::shared_ptr<PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace>> displacements,
  std::shared_ptr<PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace>> internalVirtualWork
)
{
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

template<typename Term>
void HyperelasticitySolver<Term>::
solveForDisplacements(
  std::shared_ptr<PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace>> externalVirtualWork,
  std::shared_ptr<PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace>> displacements
)
{
  // solve ∂W_int - ∂W_ext = 0 with J = 1

  // assign new value for ∂W_ext, save previous vector
  Vec rhs = externalVirtualWork->valuesGlobal();
  PetscErrorCode ierr;
  ierr = VecSwap(rhs, externalVirtualWork_); CHKERRV(ierr);

  Vec result = displacements->valuesGlobal();
  ierr = VecSwap(result, solverVariableSolution_); CHKERRV(ierr);

  nonlinearSolve();

  // restore value of ∂W_ext and solution
  ierr = VecSwap(rhs, externalVirtualWork_); CHKERRV(ierr);
  ierr = VecSwap(result, solverVariableSolution_); CHKERRV(ierr);
}

}  // namespace SpatialDiscretization
