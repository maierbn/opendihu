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
  parseSolverTypes();

  // set solver type
  ierr = KSPSetType(*ksp_, kspType_); CHKERRV(ierr);

  // set options from command line, this overrides the python config
  ierr = KSPSetFromOptions(*ksp_); CHKERRV(ierr);

  // extract preconditioner context
  PC pc;
  ierr = KSPGetPC(*ksp_, &pc); CHKERRV(ierr);

  // set type of preconditioner
  ierr = PCSetType(pc, pcType_); CHKERRV(ierr);

  // for multigrid set number of levels to 5
  if (pcType_ == std::string(PCGAMG))
  {
    int nLevels = 5;
    ierr = PCMGSetLevels(pc, nLevels, NULL); CHKERRV(ierr);
  }

  // set options from command line, this overrides the python config
  ierr = PCSetFromOptions(pc); CHKERRV(ierr);

  //                                    relative tol,      absolute tol,  diverg tol.,   max_iterations
  ierr = KSPSetTolerances (*ksp_, relativeTolerance_, PETSC_DEFAULT, PETSC_DEFAULT, maxIterations_); CHKERRV(ierr);

  // prepare log keys to log number of iterations and residual norm
  std::stringstream nIterationsLogKey;
  nIterationsLogKey << "nIterations_" << name_;
  nIterationsLogKey_ = nIterationsLogKey.str();
  
  std::stringstream residualNormLogKey;
  residualNormLogKey << "residualNorm_" << name_;
  residualNormLogKey_ = residualNormLogKey.str();
  
  std::stringstream nIterationsTotalLogKey;
  nIterationsTotalLogKey << "nIterationsTotal_" << name_;
  nIterationsTotalLogKey_ = nIterationsTotalLogKey.str();
}

void Linear::parseSolverTypes()
{
  // parse solver type
  std::string solverType = this->specificSettings_.getOptionString("solverType", "gmres");

  // parse preconditioner type
  std::string preconditionerType = this->specificSettings_.getOptionString("preconditionerType", "none");

  // all pc types: https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html
  pcType_ = PCNONE;
  if (preconditionerType == "jacobi")
  {
    pcType_ = PCJACOBI;
  }
  else if (preconditionerType == "sor")
  {
    pcType_ = PCSOR;
  }
  else if (preconditionerType == "lu")
  {
    pcType_ = PCLU;
  }
  else if (preconditionerType == "ilu")
  {
    pcType_ = PCILU;
  }
  else if (preconditionerType == "gamg")
  {
    pcType_ = PCGAMG;
  }

  // all ksp types: https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html#KSPType
  kspType_ = KSPGMRES;
  if (solverType == "richardson")
  {
    kspType_ = KSPRICHARDSON ;
  }
  else if (solverType == "chebyshev")
  {
    kspType_ = KSPCHEBYSHEV;
  }
  else if (solverType == "cg")
  {
    kspType_ = KSPCG;
  }
  else if (solverType == "preonly")
  {
    kspType_ = KSPPREONLY;
  }
  else if (solverType == "lu")
  {
    kspType_ = KSPPREONLY;
    pcType_ = PCLU;
  }
  else if (solverType == "cholesky")
  {
    kspType_ = KSPPREONLY;
    pcType_ = PCCHOLESKY;
  }
  else if (solverType == "gamg")
  {
    kspType_ = KSPPREONLY;
    pcType_ = PCGAMG;
  }
  else if (solverType == "jacobi")
  {
    kspType_ = KSPPREONLY;
    pcType_ = PCJACOBI;
  }
  else if (solverType == "sor")
  {
    kspType_ = KSPPREONLY;
    pcType_ = PCSOR;
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

void Linear::solve(Vec rightHandSide, Vec solution, std::string message)
{
  PetscErrorCode ierr;

  // solve the system
  ierr = KSPSolve(*ksp_, rightHandSide, solution); CHKERRV(ierr);

  // determine meta data
  int numberOfIterations = 0;
  PetscReal residualNorm = 0.0;
  int nDofsGlobal = 0;

  ierr = KSPGetIterationNumber(*ksp_, &numberOfIterations); CHKERRV(ierr);
  ierr = KSPGetResidualNorm(*ksp_, &residualNorm); CHKERRV(ierr);

  KSPConvergedReason convergedReason;
  ierr = KSPGetConvergedReason(*ksp_, &convergedReason); CHKERRV(ierr);
  ierr = VecGetSize(rightHandSide, &nDofsGlobal); CHKERRV(ierr);

  if (kspType_ == KSPPREONLY && (pcType_ == PCLU || pcType_ == PCILU))
  {
    if (!residual_)
    {
      temporaryVectorLeft_ = std::make_shared<Vec>();     ///< temporary vector for computation of residual for direct solvers
      temporaryVectorRight_ = std::make_shared<Vec>();    ///< temporary vector for computation of residual for direct solvers
      residual_ = std::make_shared<Vec>();    ///< residual vector for direct solvers

      Mat systemMatrix;
      ierr = KSPGetOperators(*ksp_, &systemMatrix, NULL); CHKERRV(ierr);
      ierr = MatCreateVecs(systemMatrix, &(*temporaryVectorRight_), &(*temporaryVectorLeft_)); CHKERRV(ierr);

      LOG(DEBUG) << "create temporary vectors";
    }
    LOG(DEBUG) << "compute residual";

    // compute residual
    ierr = KSPBuildResidual(*ksp_, *temporaryVectorLeft_, *temporaryVectorRight_, &(*residual_)); CHKERRV(ierr);

    LOG(INFO) << "r: " << PetscUtility::getStringVector(*residual_);

    // compute norm of residual
    ierr = VecNorm(*residual_, NORM_2, &residualNorm); CHKERRV(ierr);
  }

  if (message != "")
  {
    LOG(INFO) << message << " in " << numberOfIterations << " iterations, " << nDofsGlobal << " dofs, residual norm " << residualNorm
      << ": " << PetscUtility::getStringLinearConvergedReason(convergedReason);
  }

  // store parameter values to be logged
  Control::PerformanceMeasurement::setParameter(nIterationsLogKey_, numberOfIterations);
  Control::PerformanceMeasurement::setParameter(residualNormLogKey_, residualNorm);
  Control::PerformanceMeasurement::countNumber(nIterationsTotalLogKey_, numberOfIterations);
}

}   //namespace
