#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>

#include "control/python_config/python_config.h"

namespace Solver
{

/**
 * A solver object contains all parameters for PETSc and the PETSc solver object.
 */
class Solver
{
public:
  //! construct solver from python settings
  Solver(PythonConfig specificSettings, std::string name);

  //! destructor
  virtual ~Solver() {}

  //! determine if the own python config object is the same as config
  bool configEquals(PythonConfig config);

  //! get the name of the solver
  std::string name();

protected:

  PythonConfig specificSettings_;   //< the python config dict
  std::string name_;           //< the name of the solver
  std::string durationLogKey_;         //< key for logging of the duration of solve
};

}  // namespace
