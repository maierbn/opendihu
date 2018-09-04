#include "solver/linear.h"

#include "utility/python_utility.h"

namespace Solver
{

Linear::Linear(PyObject *specificSettings, MPI_Comm mpiCommunicator) :
  Solver(specificSettings)
{
  if (VLOG_IS_ON(1))
  {
    int size;
    PetscErrorCode ierr;
    ierr = MPI_Comm_size(mpiCommunicator, &size); CHKERRV(ierr);
    VLOG(1) << "Create linear solver on " << size << (size == 1? " rank." : " ranks.");
  }

  // parse options
  relativeTolerance_ = PythonUtility::getOptionDouble(specificSettings, "relativeTolerance", 1e-5, PythonUtility::Positive);
  maxIterations_ = PythonUtility::getOptionDouble(specificSettings, "maxIterations", 10000, PythonUtility::Positive);
  
  // set up KSP object
  //KSP *ksp;
  ksp_ = std::make_shared<KSP>();
  PetscErrorCode ierr = KSPCreate (mpiCommunicator, ksp_.get()); CHKERRV(ierr);

  // set options from command line as specified by PETSc
  KSPSetFromOptions(*ksp_);


  // extract preconditioner context
  PC pc;
  ierr = KSPGetPC (*ksp_, &pc); CHKERRV(ierr);

  // set preconditioner type
  ierr = PCSetType (pc, PCNONE); CHKERRV(ierr);

  // set solver type
  ierr = KSPSetType(*ksp_, KSPGMRES); CHKERRV(ierr);

  //                                    relative tol,      absolute tol,  diverg tol.,   max_iterations
  ierr = KSPSetTolerances (*ksp_, relativeTolerance_, PETSC_DEFAULT, PETSC_DEFAULT, maxIterations_); CHKERRV(ierr);
}

std::shared_ptr<KSP> Linear::ksp()
{
  return ksp_;
}

};   //namespace
