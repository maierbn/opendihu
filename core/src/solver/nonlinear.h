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

  //! do nothing for the nonlinear solver, the initialization is already done in the constructor
  void initialize();

protected:
  std::shared_ptr<SNES> snes_;   //< the PETSc SNES (scalable nonlinear equations solvers) object

  double snesAbsoluteTolerance_;         //< absolute solver tolerance
  double snesRelativeTolerance_;         //< relative solver tolerance
  long int snesMaxIterations_;           //< maximum number of iterations
  long int snesMaxFunctionEvaluations_;  //< maximum number of function evaluations
  int snesRebuildJacobianFrequency_;     //< how often the jacobian will be rebuild, -1 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is computed within a single nonlinear solve, 2 means every second time the Jacobian is built etc. -2 means rebuild at next chance but then never again 
  std::string snesLineSearchType_;       //< linesearch type of the snes object (SNESLineSearchType)
};

}  // namespace
