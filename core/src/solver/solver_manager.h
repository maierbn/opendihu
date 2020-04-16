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
  Manager(PythonConfig specificSettings);

  //! return previously created solver or create on the fly
  template<typename SolverType>
  std::shared_ptr<SolverType> solver(PythonConfig settings, MPI_Comm mpiCommunicator);

  //! return previously created solver or create on the fly, providing local node positions to the preconditioner
  template<typename SolverType>
  std::shared_ptr<SolverType> solver(PythonConfig settings, MPI_Comm mpiCommunicator, const std::vector<Vec3> &nodePositionsLocal);

  //! check if a solver with the given name and mpiCommunicator is stored
  bool hasSolver(std::string solverName, MPI_Comm mpiCommunicator);

  //! delete the solver identified by solverName
  void deleteSolver(std::string solverName, MPI_Comm mpiCommunicator);

  //! delete the solver identified by solverName, for all communcatiors
  void deleteSolver(std::string solverName);

private:
  //! store settings for all solvers that are specified in specificSettings_
  void storePreconfiguredSolvers();

  PythonConfig specificSettings_;    ///< python object containing the value of the python config dict with corresponding key
  int numberAnonymousSolvers_;     ///< how many inline solvers without a given name in the python config are contained in solvers_. These have a key "anonymous<no>"

  std::map<std::string, PythonConfig> solverConfiguration_;         ///< the python dicts for the solvers that were defined under "Solvers"
  std::map<MPI_Comm, std::map<std::string, std::shared_ptr<Solver>>> solvers_;    ///< for the mpi communcator the solvers with their string key
};

} // namespace
#include "solver/solver_manager.tpp"
