#include "specialized_solver/hyperelasticity/petsc_callbacks.h"

#include <Python.h>  // has to be the first included header
#include <easylogging++.h>

#include "utility/petsc_utility.h"

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

  // compute the lhs which is the virtual work and the incompressibility constraint
  object->evaluateNonlinearFunction(x, f);

  // set rows in f to x - x0 for which dirichlet BC in x is given
  object->applyDirichletBoundaryConditionsInNonlinearFunction(x, f);

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
  VLOG(1) << "pointer value jac: " << jac;
  VLOG(1) << "pointer value b:   " << b;

  object->evaluateAnalyticJacobian(x, jac);

  // zero rows and columns for which Dirichlet BC is set, set diagonal to 1
  object->applyDirichletBoundaryConditionsInJacobian(x, jac);

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
  VLOG(1) << "pointer value jac: " << jac;
  VLOG(1) << "pointer value b:   " << b;

  // compute jacobian by finite differences
  SNESComputeJacobianDefault(snes, x, jac, b, context);

  // zero rows and columns for which Dirichlet BC is set, set diagonal to 1
  object->applyDirichletBoundaryConditionsInJacobian(x, jac);

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
  VLOG(1) << "pointer value jac: " << jac;
  VLOG(1) << "pointer value b:   " << b;

  // compute the finite differences jacobian in the main jacobian slot jac
  SNESComputeJacobianDefault(snes, x, jac, jac, context);

  // zero rows and columns for which Dirichlet BC is set
  object->applyDirichletBoundaryConditionsInJacobian(x, jac);

  // compute the tangent stiffness matrix, stored in the preconditioner slot b
  object->evaluateAnalyticJacobian(x, b);

  // zero rows and columns for which Dirichlet BC is set
  object->applyDirichletBoundaryConditionsInJacobian(x, b);

  //LOG(INFO) << "terminate in jacobianFunctionCombined";
  //exit(0);
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
  //T* object = static_cast<T*>(mctx);
  LOG(DEBUG) << "  Nonlinear solver: iteration " << its << ", residual norm " << norm;

  // if log file was given, write residual norm to log file
  if (mctx != NULL)
  {
    std::ofstream *logFile = (std::ofstream *)(mctx);
    *logFile << its << ";" << norm << std::endl;
  }

  return 0;
}
