#pragma once

#include <Python.h>  // has to be the first included header
#include "solver/solver.h"

#include <memory>
#include <petscsnes.h>

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
  
  //! return the SNES object that is used for solving
  std::shared_ptr<SNES> snes();
  
protected:
  std::shared_ptr<SNES> snes_;   ///< the PETSc SNES (scalable nonlinear equations solvers) object
  double relativeTolerance_;    ///< relative solver tolerance
};

}  // namespace
