#include "solver/solver_manager.h"

namespace Solver
{
  
Manager::Manager(PyObject *specificSettings) : 
  specificSettings_(specificSettings), numberAnonymousSolvers_(0)
{
  storePreconfiguredSolvers();
}
  
bool Manager::hasSolver(std::string solverName)
{
  return solvers_.find(solverName) != solvers_.end();
}
  
void Manager::storePreconfiguredSolvers()
{
  LOG(TRACE) << "SolverManager::storePreconfiguredSolvers";
  if (specificSettings_)
  {
    if (PythonUtility::hasKey(specificSettings_, "Solvers"))
    {
      
      std::string keyString("Solvers");
      std::pair<std::string, PyObject *> dictItem 
        = PythonUtility::getOptionDictBegin<std::string, PyObject *>(specificSettings_, keyString);
      
      for (; !PythonUtility::getOptionDictEnd(specificSettings_, keyString); 
          PythonUtility::getOptionDictNext<std::string, PyObject *>(specificSettings_, keyString, dictItem))
      {
        std::string key = dictItem.first;
        PyObject *value = dictItem.second;
            
        if (value == NULL)
        {
          LOG(WARNING) << "Could not extract dict for solver \""<<key<<"\".";
        }
        else if(!PyDict_Check(value))
        {
          LOG(WARNING) << "Value for solver with name \""<<key<<"\" should be a dict.";
        }
        else
        {
          LOG(DEBUG) << "store solver configuration with key \""<<key<<"\".";
          solverConfiguration_[key] = value;
        }
      }
    }
    else
    {
      LOG(INFO) << "You have specified the solver in-line and not under the extra key \"Solvers\". You could do so, "
        " by defining \"Solvers\": {\"<your custom solver name>\": {<your solver parameters>}} at the beginning of the "
        " config and \"solverName\": \"<your custom solver name>\" where you currently have specified the solver parameters. "
        " This is required if you want to use the same solver for multiple objects.";
    }
  }
}
  
};   // namespace