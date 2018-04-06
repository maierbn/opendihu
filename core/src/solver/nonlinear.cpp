#include "solver/nonlinear.h"

#include "utility/python_utility.h"

namespace Solver
{
  
Nonlinear::Nonlinear(PyObject *specificSettings) : Solver(specificSettings)
{
  // parse options
  relativeTolerance_ = PythonUtility::getOptionDouble(specificSettings, "relativeTolerance", 1e-5, PythonUtility::Positive);

  // set up SNES object
  snes_ = std::make_shared<SNES>();
  PetscErrorCode ierr = SNESCreate (PETSC_COMM_WORLD, snes_.get()); CHKERRV(ierr);
  
  // set options from command line as specified by PETSc
  SNESSetFromOptions(*snes_);
  
  // extract linear solver context
  KSP ksp;
  ierr = SNESGetKSP (*snes_, &ksp); CHKERRV(ierr);
  
  // set solver type
  ierr = KSPSetType(ksp, KSPGMRES); CHKERRV(ierr);
  
  //                             relative tol,       absolute tol,  diverg tol.,   max_iterations
  ierr = KSPSetTolerances (ksp, relativeTolerance_, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRV(ierr);
}
  
std::shared_ptr<SNES> Nonlinear::snes()
{
  return snes_;
}

};   //namespace