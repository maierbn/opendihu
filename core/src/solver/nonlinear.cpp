#include "solver/nonlinear.h"

#include "utility/python_utility.h"

namespace Solver
{

Nonlinear::Nonlinear(PythonConfig specificSettings, MPI_Comm mpiCommunicator, std::string name) :
  Linear(specificSettings, mpiCommunicator, name)
{
  // set up SNES object
  snes_ = std::make_shared<SNES>();
  PetscErrorCode ierr = SNESCreate (mpiCommunicator, snes_.get()); CHKERRV(ierr);

  // set options from command line as specified by PETSc
  SNESSetFromOptions(*snes_);

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
