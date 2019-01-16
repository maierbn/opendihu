#include "solver/linear.h"

#include "utility/python_utility.h"
#include "control/performance_measurement.h"

namespace Solver
{

Linear::Linear(PythonConfig specificSettings, MPI_Comm mpiCommunicator, std::string name) :
  Solver(specificSettings, name)
{
  if (VLOG_IS_ON(1))
  {
    int size;
    PetscErrorCode ierr;
    ierr = MPI_Comm_size(mpiCommunicator, &size); CHKERRV(ierr);
    VLOG(1) << "Create linear solver on " << size << (size == 1? " rank." : " ranks.");
  }

  // parse options
  relativeTolerance_ = this->specificSettings_.getOptionDouble("relativeTolerance", 1e-5, PythonUtility::Positive);
  maxIterations_ = this->specificSettings_.getOptionDouble("maxIterations", 10000, PythonUtility::Positive);
  
  // set up KSP object
  //KSP *ksp;
  ksp_ = std::make_shared<KSP>();
  PetscErrorCode ierr = KSPCreate (mpiCommunicator, ksp_.get()); CHKERRV(ierr);

  // parse the solver and preconditioner types from settings
  KSPType kspType;
  PCType pcType;
  parseSolverTypes(kspType, pcType);

  // set solver type
  ierr = KSPSetType(*ksp_, kspType); CHKERRV(ierr);

  // set options from command line, this overrides the python config
  ierr = KSPSetFromOptions(*ksp_); CHKERRV(ierr);

  // extract preconditioner context
  PC pc;
  ierr = KSPGetPC(*ksp_, &pc); CHKERRV(ierr);

  // set type of preconditioner
  ierr = PCSetType(pc, pcType); CHKERRV(ierr);

  // for multigrid set number of levels to 5
  if (pcType == std::string(PCGAMG))
  {
    int nLevels = 5;
    ierr = PCMGSetLevels(pc, nLevels, NULL); CHKERRV(ierr);
  }

  // set options from command line, this overrides the python config
  ierr = PCSetFromOptions(pc); CHKERRV(ierr);

  //                                    relative tol,      absolute tol,  diverg tol.,   max_iterations
  ierr = KSPSetTolerances (*ksp_, relativeTolerance_, PETSC_DEFAULT, PETSC_DEFAULT, maxIterations_); CHKERRV(ierr);
}

void Linear::parseSolverTypes(KSPType &kspType, PCType &pcType)
{
  // parse solver type
  std::string solverType = this->specificSettings_.getOptionString("solverType", "gmres");

  // parse preconditioner type
  std::string preconditionerType = this->specificSettings_.getOptionString("preconditionerType", "none");

  // all pc types: https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html
  pcType = PCNONE;
  if (preconditionerType == "jacobi")
  {
    pcType = PCJACOBI;
  }
  else if (preconditionerType == "sor")
  {
    pcType = PCSOR;
  }
  else if (preconditionerType == "lu")
  {
    pcType = PCLU;
  }
  else if (preconditionerType == "ilu")
  {
    pcType = PCILU;
  }
  else if (preconditionerType == "gamg")
  {
    pcType = PCGAMG;
  }

  // all ksp types: https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html#KSPType
  kspType = KSPGMRES;
  if (solverType == "richardson")
  {
    kspType = KSPRICHARDSON ;
  }
  else if (solverType == "chebyshev")
  {
    kspType = KSPCHEBYSHEV;
  }
  else if (solverType == "cg")
  {
    kspType = KSPCG;
  }
  else if (solverType == "preonly")
  {
    kspType = KSPPREONLY;
  }
  else if (solverType == "lu")
  {
    kspType = KSPPREONLY;
    pcType = PCLU;
  }
  else if (solverType == "cholesky")
  {
    kspType = KSPPREONLY;
    pcType = PCCHOLESKY;
  }
  else if (solverType == "gamg")
  {
    kspType = KSPPREONLY;
    pcType = PCGAMG;
  }

  std::stringstream optionKey;
  optionKey << this->name_ << "_solverType";
  Control::PerformanceMeasurement::setParameter(optionKey.str(), solverType);

  optionKey.str("");
  optionKey << this->name_ << "_preconditionerType";
  Control::PerformanceMeasurement::setParameter(optionKey.str(), preconditionerType);
}

std::shared_ptr<KSP> Linear::ksp()
{
  return ksp_;
}

LinearMG::LinearMG(PythonConfig specificSettings, MPI_Comm mpiCommunicator, std::string name) :
  Solver(specificSettings, name)
{
  if (VLOG_IS_ON(1))
  {
    int size;
    PetscErrorCode ierr;
    ierr = MPI_Comm_size(mpiCommunicator, &size); CHKERRV(ierr);
    VLOG(1) << "Create linear solver on " << size << (size == 1? " rank." : " ranks.");
  }

  // parse options
  relativeTolerance_ = this->specificSettings_.getOptionDouble("relativeTolerance", 1e-5, PythonUtility::Positive);
  maxIterations_ = this->specificSettings_.getOptionDouble("maxIterations", 10000, PythonUtility::Positive);
  
  // set up KSP object
  //KSP *ksp;
  ksp_ = std::make_shared<KSP>();
  PetscErrorCode ierr = KSPCreate (mpiCommunicator, ksp_.get()); CHKERRV(ierr);
  
  ksp2_ = std::make_shared<KSP>();
  PetscErrorCode ierr2 = KSPCreate (mpiCommunicator, ksp2_.get()); CHKERRV(ierr2);
  
  
  // parse the solver and preconditioner types from settings
  KSPType kspType;
  PCType pcType;
  parseSolverTypes(kspType, pcType);

  // set solver type
  ierr = KSPSetType(*ksp_, kspType); CHKERRV(ierr);
  ierr2 = KSPSetType(*ksp2_, kspType); CHKERRV(ierr2);
  
  // set options from command line, this overrides the python config
  ierr = KSPSetFromOptions(*ksp_); CHKERRV(ierr);
  ierr2 = KSPSetFromOptions(*ksp2_); CHKERRV(ierr2);
  
  // extract preconditioner context
  PC pc;
  ierr = KSPGetPC(*ksp_, &pc); CHKERRV(ierr);
  ierr2 = KSPGetPC(*ksp2_, &pc); CHKERRV(ierr2);

  // set type of preconditioner
  ierr = PCSetType(pc, pcType); CHKERRV(ierr);
  ierr2 = PCSetType(pc, pcType); CHKERRV(ierr2);
  
  // for multigrid set number of levels to 5
  //TODO ???
  if (pcType == std::string(PCGAMG))
  {
    int nLevels = 5;
    ierr = PCMGSetLevels(pc, nLevels, NULL); CHKERRV(ierr);
  }

  // set options from command line, this overrides the python config
  ierr = PCSetFromOptions(pc); CHKERRV(ierr);
  ierr2 = PCSetFromOptions(pc); CHKERRV(ierr2);

  //                                    relative tol,      absolute tol,  diverg tol.,   max_iterations
  ierr = KSPSetTolerances (*ksp_, relativeTolerance_, PETSC_DEFAULT, PETSC_DEFAULT, maxIterations_); CHKERRV(ierr);
  ierr2 = KSPSetTolerances (*ksp2_, relativeTolerance_, PETSC_DEFAULT, PETSC_DEFAULT, maxIterations_); CHKERRV(ierr2);
}

void LinearMG::parseSolverTypes(KSPType &kspType, PCType &pcType)
{
  // parse solver type
  std::string solverType = this->specificSettings_.getOptionString("solverType", "gmres");

  // parse preconditioner type
  std::string preconditionerType = this->specificSettings_.getOptionString("preconditionerType", "none");

  // all pc types: https://www.mcs.anl.goierr = PCSetType(pc, pcType); CHKERRV(ierr);v/petsc/petsc-current/docs/manualpages/PC/PCType.html
  pcType = PCNONE;
  if (preconditionerType == "jacobi")
  {
    pcType = PCJACOBI;
  }
  else if (preconditionerType == "sor")
  {
    pcType = PCSOR;
  }
  else if (preconditionerType == "lu")
  {
    pcType = PCLU;
  }
  else if (preconditionerType == "ilu")
  {
    pcType = PCILU;
  }
  else if (preconditionerType == "gamg")
  {
    pcType = PCGAMG;
  }

  // all ksp types: https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html#KSPType
  kspType = KSPGMRES;
  if (solverType == "richardson")
  {
    kspType = KSPRICHARDSON ;
  }
  else if (solverType == "chebyshev")
  {
    kspType = KSPCHEBYSHEV;
  }
  else if (solverType == "cg")
  {
    kspType = KSPCG;
  }
  else if (solverType == "preonly")
  {
    kspType = KSPPREONLY;
  }
  else if (solverType == "lu")
  {
    kspType = KSPPREONLY;
    pcType = PCLU;
  }
  else if (solverType == "cholesky")
  {
    kspType = KSPPREONLY;
    pcType = PCCHOLESKY;
  }
  else if (solverType == "gamg")
  {
    kspType = KSPPREONLY;
    pcType = PCGAMG;
  }

  std::stringstream optionKey;
  optionKey << this->name_ << "_solverType";
  Control::PerformanceMeasurement::setParameter(optionKey.str(), solverType);

  optionKey.str("");
  optionKey << this->name_ << "_preconditionerType";
  Control::PerformanceMeasurement::setParameter(optionKey.str(), preconditionerType);
}
std::shared_ptr<KSP> LinearMG::ksp()
{
  return ksp_;
}
std::shared_ptr<KSP> LinearMG::ksp2()
{
  return ksp2_;
}
std::shared_ptr<int> LinearMG::numCycles()
{
  return numCycles_;
}
};   //namespace
