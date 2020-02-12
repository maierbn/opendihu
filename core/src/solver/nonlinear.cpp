#include "solver/nonlinear.h"

#include "utility/python_utility.h"

namespace Solver
{

Nonlinear::Nonlinear(PythonConfig specificSettings, MPI_Comm mpiCommunicator, std::string name) :
  Linear(specificSettings, mpiCommunicator, name)
{
  snesAbsoluteTolerance_ = this->specificSettings_.getOptionDouble("snesAbsoluteTolerance", 1e-10, PythonUtility::Positive);
  snesRelativeTolerance_ = this->specificSettings_.getOptionDouble("snesRelativeTolerance", 1e-10, PythonUtility::Positive);
  snesMaxIterations_ = this->specificSettings_.getOptionDouble("snesMaxIterations", 50, PythonUtility::Positive);
  snesMaxFunctionEvaluations_ = this->specificSettings_.getOptionDouble("snesMaxFunctionEvaluations", 1000, PythonUtility::Positive);
  snesLineSearchType_ = this->specificSettings_.getOptionString("snesLineSearchType", "l2");

  // assert that snesLineSearchType_ has a valid type: https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESLineSearchType.html#SNESLineSearchType
  if (snesLineSearchType_ != "bt" && snesLineSearchType_ != "nleqerr" && snesLineSearchType_ != "basic" && snesLineSearchType_ != "l2" && snesLineSearchType_ != "cp" && snesLineSearchType_ != "shell" && snesLineSearchType_ != "ncglinear")
  {
    LOG(ERROR) << this->specificSettings_ << "[\"snesLineSearchType\"] has invalid value \"" << snesLineSearchType_ << "\". "
      << "Allowed values are \"bt\" \"nleqerr\" \"basic\" \"l2\" \"cp\" \"ncglinear\". Now using default value \"l2\".";
    snesLineSearchType_ = "l2";
  }

  // set up SNES object
  snes_ = std::make_shared<SNES>();
  PetscErrorCode ierr = SNESCreate (mpiCommunicator, snes_.get()); CHKERRV(ierr);

  // PetscErrorCode  SNESSetTolerances(SNES snes,PetscReal abstol,PetscReal rtol,PetscReal stol,PetscInt maxit,PetscInt maxf)
  ierr = SNESSetTolerances(*snes_, snesAbsoluteTolerance_, snesRelativeTolerance_, PETSC_DEFAULT, snesMaxIterations_, snesMaxFunctionEvaluations_); CHKERRV(ierr);

  // set options from command line as specified by PETSc
  ierr = SNESSetFromOptions(*snes_); CHKERRV(ierr);

  // extract linesearch context
  SNESLineSearch linesearch;
  ierr = SNESGetLineSearch(*snes_, &linesearch); CHKERRV(ierr);

  //! set line search type from settings
  ierr = SNESLineSearchSetType(linesearch, (SNESLineSearchType)snesLineSearchType_.c_str()); CHKERRV(ierr);

  //! override line search settings from command line
  ierr = SNESLineSearchSetFromOptions(linesearch); CHKERRV(ierr);

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
