#include "solver/solver_manager.h"

#include <Python.h> // this header has to be included first

#include <omp.h>
#include <memory>
#include <iostream>
#include "easylogging++.h"
#include "utility/python_utility.h"
#include "utility/string_utility.h"

namespace Solver
{

//! return previously created solver or create on the fly
template<typename SolverType>
std::shared_ptr<SolverType> Manager::solver(PythonConfig settings, MPI_Comm mpiCommunicator, std::string solverNameKey)
{
  LOG(TRACE) << "Manager::solver";

  // if there is not yet an entry for the mpi communicator, create an empty one
  if (solvers_.find(mpiCommunicator) == solvers_.end())
  {
    solvers_[mpiCommunicator] = std::map<std::string, std::shared_ptr<Solver>>();
  }

  // if solver has already been created earlier
  if (settings.hasKey(solverNameKey))
  {
    std::string solverName = settings.getOptionString(solverNameKey, "");
    
    if (hasSolver(solverName, mpiCommunicator))
    {
      LOG(DEBUG) << "Solver with solverName \"" << solverName << "\" requested and found, type is "
        << StringUtility::demangle(typeid(solvers_[mpiCommunicator][solverName]).name());

      return std::static_pointer_cast<SolverType>(solvers_[mpiCommunicator][solverName]);
    }
    else if (solverConfiguration_.find(solverName) != solverConfiguration_.end())
    {
      // solver was preconfigured, do nothing specific here, created standard solver
      LOG(DEBUG) << "Solver configuration for \"" << solverName << "\" requested and found, create solver. "
        << "Type is " << StringUtility::demangle(typeid(SolverType).name()) << ".";

      // create new solver object
      std::shared_ptr<SolverType> solver = std::make_shared<SolverType>(solverConfiguration_.at(solverName), mpiCommunicator, solverName);
      
      // initialize solver
      solver->initialize();

      solvers_[mpiCommunicator][solverName] = solver;

      LOG(DEBUG) << "Stored under key \"" << solverName << "\".";
      return std::static_pointer_cast<SolverType>(solvers_[mpiCommunicator][solverName]);
    }
    else
    {
      LOG(ERROR) << "Config contains reference to solver with solverName \"" << solverName << "\" but no such solver was defined.";
    }
  }
  else
  {
    VLOG(1) << "Config does not contain solverName.";
  }

  // check if there is a matching solver already stored
  // loop over all stored solvers
  for (auto &solver: solvers_[mpiCommunicator])
  {
    // check if type matches
    if (std::dynamic_pointer_cast<SolverType>(solver.second))
    {
      // check if config is the  same
      if (solver.second->configEquals(settings))
      {
        LOG(DEBUG) << "Solver \"" << solver.first << "\" matches settings.";
        VLOG(1) << "Solver \"" << solver.first << "\" matches settings.";
        return std::static_pointer_cast<SolverType>(solver.second);
      }
    }
  }


  // create new solver, store as anonymous object
  std::stringstream anonymousName;
  anonymousName << "anonymous" << numberAnonymousSolvers_++;
  LOG(DEBUG) << "Create new solver with type " << StringUtility::demangle(typeid(SolverType).name())
    << " and name \"" <<anonymousName.str() << "\".";
  
  // create new solver object
  std::shared_ptr<SolverType> solver = std::make_shared<SolverType>(settings, mpiCommunicator, anonymousName.str());
  
  // initialize solver
  solver->initialize();

  solvers_[mpiCommunicator][anonymousName.str()] = solver;

  return solver;
}

}  // namespace
