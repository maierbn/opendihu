#pragma once

#include <Python.h>  // has to be the first included header

#include <iostream>

namespace Solver
{

/**
 * A solver object contains all parameters for PETSc and the PETSc solver object.
 */
class Solver
{
public:
  //! construct solver from python settings
  Solver(PyObject *specificSettings, std::string name);
  virtual ~Solver() {}

  //! determine if the own python config object is the same as config
  bool configEquals(PyObject *config);
protected:

  PyObject *specificSettings_;   ///< the python config dict
  std::string name_;           ///< the name of the solver
};

}  // namespace
