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
  Linear(PyObject *specificSettings);
  
  //! return the KSP object that is used for solving
  std::shared_ptr<KSP> ksp();
  
protected:
 
  std::shared_ptr<KSP> ksp_;   ///< the PETSc KSP (Krylov subspace) object
  double relativeTolerance_;    ///< relative solver tolerance
};

}  // namespace
