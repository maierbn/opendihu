#pragma once

#include <Python.h>  // has to be the first included header
#include "solver/solver.h"

#include <petscksp.h>
#include <memory>

namespace Solver
{

/**
 * A linear solver configuration that contains all parameters for PETSc and the KSP object.
 */
class Linear : public Solver
{
public:
  //! construct solver from python settings
  Linear(PythonConfig specificSettings, MPI_Comm mpiCommunicator, std::string name);

  //! return the KSP object that is used for solving
  std::shared_ptr<KSP> ksp();

  //! perform the solve
  void solve(Vec rightHandSide, Vec solution, std::string message="");

protected:

  //! parse the solver and preconditioner type from settings
  void parseSolverTypes();

  std::shared_ptr<KSP> ksp_;   ///< the PETSc KSP (Krylov subspace) object
  double relativeTolerance_;    ///< relative solver tolerance
  long int maxIterations_;     ///< maximum number of iterations

  KSPType kspType_;    ///< the solver type
  PCType pcType_;      ///< the preconditioner type

  std::shared_ptr<Vec> temporaryVectorLeft_;     ///< temporary vector for computation of residual for direct solvers
  std::shared_ptr<Vec> temporaryVectorRight_;    ///< temporary vector for computation of residual for direct solvers
  std::shared_ptr<Vec> residual_;    ///< residual vector for direct solvers


  std::string nIterationsLogKey_;  ///< the keyword for the log with which the number of iterations will be stored
  std::string residualNormLogKey_;  ///< the keyword for the log with which the residual norm gets stored
};

}  // namespace
