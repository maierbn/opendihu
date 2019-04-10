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
  void parseSolverTypes(KSPType &kspType, PCType &pcType);

  std::shared_ptr<KSP> ksp_;   ///< the PETSc KSP (Krylov subspace) object
  double relativeTolerance_;    ///< relative solver tolerance
  long int maxIterations_;     ///< maximum number of iterations

  std::string nIterationsLogKey_;  ///< the keyword for the log with which the number of iterations will be stored
  std::string residualNormLogKey_;  ///< the keyword for the log with which the residual norm gets stored
  std::string nIterationsTotalLogKey_;  ///< the keyword for the log with which the total number of iterations gets stored
};

}  // namespace
