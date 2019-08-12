#pragma once

#include "petsc.h"

//! compute the numeric jacobian, this is an adjusted copy of the original PETSc implementation with little modifications such that it also works
//! with Nested vectors
PetscErrorCode  SNESComputeJacobianDefaultNested(SNES snes,Vec x1,Mat J,Mat B,void *ctx);
