#include "specialized_solver/solid_mechanics/hyperelasticity/petsc_callbacks.h"

#include <Python.h>  // has to be the first included header
#include <easylogging++.h>

#include "utility/petsc_utility.h"
#include "specialized_solver/solid_mechanics/hyperelasticity/compute_numeric_jacobian.h"

/**
 * Nonlinear function F that gets solved by PETSc, solve for x such that F(x) = 0
 *  Input Parameters:
 *  snes - the SNES context
 *  x    - input vector
 *  context  - optional user-defined context
 *
 *  Output Parameter:
 *  f - function vector
 */
template<typename T>
PetscErrorCode nonlinearFunction(SNES snes, Vec x, Vec f, void *context)
{
  T* object = static_cast<T*>(context);

  VLOG(1) << "in nonlinearFunction";
  VLOG(1) << "pointer value x: " << x;
  VLOG(1) << "pointer value f: " << f;

  // set prescribed values in x to x0, to also make the columns vanish in the numeric jacobian
  //object->applyDirichletBoundaryConditionsInVector(x);

  // compute the lhs which is the virtual work and the incompressibility constraint
  object->evaluateNonlinearFunction(x, f);

  // set rows in f to x - x0 for which dirichlet BC in x is given
  //object->applyDirichletBoundaryConditionsInNonlinearFunction(x_original, f);

  // compute and output function norm
  if (VLOG_IS_ON(1))
  {
    PetscReal functionNorm;
    VecNorm(f, NORM_2, &functionNorm);
    VLOG(1) << "function norm: " << functionNorm;
  }
  return 0;
}

/**
 * Evaluates Jacobian matrix
 *  Input Parameters:
 *  snes - the SNES context
 *  x    - input vector
 *  context  - optional user-defined context
 *
 *  Output Parameters:
 *  jac - Jacobian matrix
 *  b   - optionally different preconditioning matrix
 */
template<typename T>
PetscErrorCode jacobianFunctionAnalytic(SNES snes, Vec x, Mat jac, Mat b, void *context)
{
  T* object = static_cast<T*>(context);

  VLOG(1) << "in jacobianFunctionAnalytic";
  VLOG(1) << "pointer value x:   " << x;
  VLOG(1) << "pointer value jac: " << jac << " (should be analytic slot)";
  VLOG(1) << "pointer value b:   " << b << " (should be analytic slot)";

  // compute jacobian by analytic formula
  object->evaluateAnalyticJacobian(x, jac);

  // zero rows and columns for which Dirichlet BC is set, set diagonal to 1
  object->applyDirichletBoundaryConditionsInJacobian(x, jac);

  // output the jacobian matrix for debugging
  object->dumpJacobianMatrix(jac);

  VLOG(2) << "-- computed tangent stiffness matrix analytically: " << PetscUtility::getStringMatrix(jac);
  VLOG(2) << "-- non-zeros pattern: " << std::endl << PetscUtility::getStringSparsityPattern(jac);

  return 0;
}

template<typename T>
PetscErrorCode jacobianFunctionFiniteDifferences(SNES snes, Vec x, Mat jac, Mat b, void *context)
{
  T* object = static_cast<T*>(context);

  VLOG(1) << "in jacobianFunctionFiniteDifferences";
  VLOG(1) << "pointer value x:   " << x;
  VLOG(1) << "pointer value jac: " << jac << " (should be numeric slot)";
  VLOG(1) << "pointer value b:   " << b << " (should be numeric slot)";

  LOG(DEBUG) << "in jacobianFunctionFiniteDifferences, "
    << "solution: " << object->combinedVecSolution()->getString() << ", residual: " << object->combinedVecResidual()->getString();

  // compute jacobian by finite differences, in b (but this is the same pointer as jac)
  SNESComputeJacobianDefaultNested(snes, x, jac, b, context);

  // zero rows and columns for which Dirichlet BC is set, set diagonal to 1
  object->applyDirichletBoundaryConditionsInJacobian(x, jac);

  // output the jacobian matrix for debugging
  object->dumpJacobianMatrix(jac);

  VLOG(2) << "-- computed tangent stiffness matrix by finite differences: " << PetscUtility::getStringMatrix(jac);
  VLOG(2) << "-- non-zeros pattern: " << std::endl << PetscUtility::getStringSparsityPattern(jac);

  return 0;
}


template<typename T>
PetscErrorCode jacobianFunctionCombined(SNES snes, Vec x, Mat jac, Mat b, void *context)
{
  T* object = static_cast<T*>(context);

  VLOG(1) << "in jacobianFunctionCombined";
  VLOG(1) << "pointer value x:   " << x;
  VLOG(1) << "pointer value jac: " << jac << " (should be the numeric slot)";
  VLOG(1) << "pointer value b:   " << b << " (should be the analytic slot)";

  // compute the finite differences jacobian in the main jacobian slot jac
  SNESComputeJacobianDefaultNested(snes, x, jac, jac, context);

  // zero rows and columns for which Dirichlet BC is set
  object->applyDirichletBoundaryConditionsInJacobian(x, jac);

  // output the jacobian matrix for debugging
  object->dumpJacobianMatrix(jac);

  // compute the tangent stiffness matrix, stored in the preconditioner slot b
  object->evaluateAnalyticJacobian(x, b);

  // zero rows and columns for which Dirichlet BC is set
  object->applyDirichletBoundaryConditionsInJacobian(x, b);

  // output the jacobian matrix for debugging
  object->dumpJacobianMatrix(b);

  object->debug();

  //LOG_AFTER_N(2,FATAL) << "terminate in jacobianFunctionCombined";
  return 0;
}

/**
 * Monitor convergence of nonlinear solver
 *
 *  Input Parameters:
 *   snes  - the SNES context
 *   its   - iteration number
 *   norm  - 2-norm function value (may be estimated)
 *   mctx  - [optional] monitoring context
 */
template<typename T>
PetscErrorCode monitorFunction(SNES snes, PetscInt its, PetscReal norm, void *mctx)
{
  T* object = static_cast<T*>(mctx);
  object->monitorSolvingIteration(snes, its, norm);

  return 0;
}
