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
class Newton : public Solver
{
public:
  //! construct solver from python settings
  Newton(PyObject *specificSettings);
  
  //! return the SNES object that is used for solving
  std::shared_ptr<SNES> ksp();
  
protected:
  std::shared_ptr<SNES> snes_;   ///< the PETSc SNES (scalable nonlinear equations solvers) object
};

}  // namespace
