#pragma once

#include <Python.h>  // has to be the first included header
#include "solver/solver.h"

#include <memory>
#include <petscsnes.h>
#include <petscksp.h>

namespace Solver
{

/**
 * A nonlinear Newton solver configuration that contains all parameters for PETSc and the SNES object.
 */
class Nonlinear : public Solver
{
public:
  //! construct solver from python settings
  Nonlinear(PyObject *specificSettings);

  //! return the SNES object that is used for solving (nonlinear solver context)
  std::shared_ptr<SNES> snes();

  //! return the KSP object of the linear solver that is used by the nonlinear solver
  std::shared_ptr<KSP> ksp();

protected:
  std::shared_ptr<SNES> snes_;   ///< the PETSc SNES (scalable nonlinear equations solvers) object
  std::shared_ptr<KSP> ksp_;     ///< KSP object of the linear solver that is used by the nonlinear solver

  double relativeTolerance_;    ///< relative solver tolerance
};

}  // namespace
