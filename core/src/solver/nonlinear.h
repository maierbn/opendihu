#pragma once

#include <Python.h>  // has to be the first included header
#include "solver/linear.h"

#include <memory>
#include <petscsnes.h>
#include <petscksp.h>

namespace Solver
{

/**
 * A nonlinear Newton solver configuration that contains all parameters for PETSc and the SNES object.
 */
class Nonlinear : public Linear
{
public:
  //! construct solver from python settings
  Nonlinear(PythonConfig specificSettings, MPI_Comm mpiCommunicator, std::string name);

  //! return the SNES object that is used for solving (nonlinear solver context)
  std::shared_ptr<SNES> snes();

protected:
  std::shared_ptr<SNES> snes_;   ///< the PETSc SNES (scalable nonlinear equations solvers) object

  double snesAbsoluteTolerance_;         ///< absolute solver tolerance
  double snesRelativeTolerance_;         ///< relative solver tolerance
  long int snesMaxIterations_;           ///< maximum number of iterations
  long int snesMaxFunctionEvaluations_;  ///< maximum number of function evaluations
};

}  // namespace
