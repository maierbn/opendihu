#include "solver/nonlinear.h"

#include "utility/python_utility.h"

namespace Solver
{

Nonlinear::Nonlinear(PythonConfig specificSettings, MPI_Comm mpiCommunicator, std::string name) :
  Linear(specificSettings, mpiCommunicator, name)
{
  snesRelativeTolerance_ = this->specificSettings_.getOptionDouble("snesRelativeTolerance", 1e-10, PythonUtility::Positive);
  snesMaxIterations_ = this->specificSettings_.getOptionDouble("snesMaxIterations", 50, PythonUtility::Positive);
  snesMaxFunctionEvaluations_ = this->specificSettings_.getOptionDouble("snesMaxFunctionEvaluations", 1000, PythonUtility::Positive);

  // set up SNES object
  snes_ = std::make_shared<SNES>();
  PetscErrorCode ierr = SNESCreate (mpiCommunicator, snes_.get()); CHKERRV(ierr);

  // set options from command line as specified by PETSc
  ierr = SNESSetFromOptions(*snes_); CHKERRV(ierr);

  // PetscErrorCode  SNESSetTolerances(SNES snes,PetscReal abstol,PetscReal rtol,PetscReal stol,PetscInt maxit,PetscInt maxf)
  ierr = SNESSetTolerances(*snes_, PETSC_DEFAULT, snesRelativeTolerance_, PETSC_DEFAULT, snesMaxIterations_, snesMaxFunctionEvaluations_); CHKERRV(ierr);

  // extract linear solver context
  ksp_ = std::make_shared<KSP>();
  ierr = SNESGetKSP (*snes_, ksp_.get()); CHKERRV(ierr);

  this->setupKsp(*this->ksp_);
}

std::shared_ptr<SNES> Nonlinear::snes()
{
  return snes_;
}

}   //namespace
