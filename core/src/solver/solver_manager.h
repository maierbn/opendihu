#pragma once

#include <Python.h>  // has to be the first included header
#include <map>

#include "control/dihu_context.h"
#include "solver/solver.h"

namespace Solver 
{

/**
 * This class creates and stores all used PETSc solver handles. It works analoguos to mesh_manager.
 * Each solver configuration can be defined in the python config under "Solvers" with a name and other properties.
 * The name can be referenced in the config where the solver is needed, e.g. under "FiniteElements".
 * Instead of specification under "Solvers" the parameters can also be specified directly.
 */
class Manager
{
public:
  //! constructor
  Manager(PyObject *specificSettings);
  
  //! return previously created solver or create on the fly
  template<typename SolverType>
  std::shared_ptr<SolverType> solver(PyObject *settings);
  
  //! check if a solver with the given name is stored
  bool hasSolver(std::string solverName);
  
private:
  //! store settings for all solvers that are specified in specificSettings_
  void storePreconfiguredSolvers();
 
  PyObject *specificSettings_;    ///< python object containing the value of the python config dict with corresponding key
  int numberAnonymousSolvers_;     ///< how many inline solvers without a given name in the python config are contained in solvers_. These have a key "anonymous<no>"
  
  std::map<std::string, PyObject *> solverConfiguration_;         ///< the python dicts for the solvers that were defined under "Solvers"
  std::map<std::string, std::shared_ptr<Solver>> solvers_;    ///< the solvers with their string key
};

};  // namespace
#include "solver/solver_manager.tpp"
