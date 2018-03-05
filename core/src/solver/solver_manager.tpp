#include "solver/solver_manager.h"

#include "Python.h" // this header has to be included first

#include <memory>
#include <iostream>
#include "easylogging++.h"

namespace Solver
{

//! return previously created solver or create on the fly
template<typename SolverType>
std::shared_ptr<SolverType> Manager::solver(PyObject *settings)
{
  // if solver has already been created earlier
  if (PythonUtility::containsKey(settings, "solverName"))
  {
    std::string solverName = PythonUtility::getOptionString(settings, "solverName", "");
    if (hasSolver(solverName))
    {
      VLOG(1) << "Solver with solverName \""<<solverName<<"\" requested and found, type is "<<typeid(solvers_[solverName]).name();
      return std::static_pointer_cast<SolverType>(solvers_[solverName]);
    }
    else if(solverConfiguration_.find(solverName) != solverConfiguration_.end())
    {
      // solver was preconfigured, do nothing specific here, created standard solver with 1 node
      LOG(DEBUG) << "Solver configuration for \""<<solverName<<"\" found and requested, will be created now. "
        << "Type is "<< typeid(SolverType).name()<<".";
      std::shared_ptr<SolverType> solver = std::make_shared<SolverType>(solverConfiguration_[solverName]);
      solvers_[solverName] = solver;
      LOG(DEBUG) << "Stored under key \""<<solverName<<"\".";
      return std::static_pointer_cast<SolverType>(solvers_[solverName]);
    }
    else
    {
      LOG(ERROR) << "Config contains reference to solver with solverName \""<<solverName<<"\" but no such solver was defined.";      
    }
  }
  else
  {
    LOG(DEBUG) << "Config does not contain solverName.";
  }
  
  // check if there is a matching solver already stored
  // loop over all stored solvers 
  for(auto &solver: solvers_)
  {
    // check if type matches
    if (std::dynamic_pointer_cast<SolverType>(solver.second))
    {
      // check if config is the  same
      if (solver.second->configEquals(settings))
      {
        VLOG(1) << "Solver \"" << solver.first << "\" matches settings.";
        return std::static_pointer_cast<SolverType>(solver.second);
      }
    }
  }
  
  
  // create new solver, store as anonymous object
  std::stringstream anonymousName;
  anonymousName << "anonymous" << numberAnonymousSolvers_++;
  LOG(DEBUG) << "Create new solver with type "<<typeid(SolverType).name()<<" and name \""<<anonymousName.str()<<"\".";
  std::shared_ptr<SolverType> solver = std::make_shared<SolverType>(settings);

  solvers_[anonymousName.str()] = solver;
  
  return solver;
}
  
};   // namespace 