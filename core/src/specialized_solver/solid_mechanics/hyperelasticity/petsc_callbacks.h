#pragma once

#include <petsc/petsc.h>

/**
 * Callback functions for PETSc and the SNES solver
 */

/**
 * Nonlinear function F that gets solved by PETSc, solve for u such that F(u) = 0
 *  Input Parameters:
 *  snes - the SNES context
 *  u    - input vector
 *  context  - optional user-defined context
 *
 *  Output Parameter:
 *  f - function vector
 */
template<typename T>
PetscErrorCode nonlinearFunction(SNES snes, Vec u, Vec f, void *context);

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
PetscErrorCode jacobianFunctionAnalytic(SNES snes, Vec x, Mat jac, Mat b, void *context);


template<typename T>
PetscErrorCode jacobianFunctionFiniteDifferences(SNES snes, Vec x, Mat jac, Mat b, void *context);


template<typename T>
PetscErrorCode jacobianFunctionCombined(SNES snes, Vec x, Mat jac, Mat b, void *context);

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
PetscErrorCode monitorFunction(SNES snes, PetscInt its, PetscReal norm, void *mctx);

#include "specialized_solver/solid_mechanics/hyperelasticity/petsc_callbacks.tpp"
