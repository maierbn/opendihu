#include "solver/solver_manager.h"

namespace Solver
{

Manager::Manager(PythonConfig specificSettings) :
  specificSettings_(specificSettings), numberAnonymousSolvers_(0)
{
  storePreconfiguredSolvers();
}

bool Manager::hasSolver(std::string solverName, MPI_Comm mpiCommunicator)
{
  // check if solvers for the given mpiCommuncator exist
  if (solvers_.find(mpiCommunicator) == solvers_.end())
  {
    return false;
  }

  return solvers_[mpiCommunicator].find(solverName) != solvers_[mpiCommunicator].end();
}

void Manager::storePreconfiguredSolvers()
{
  LOG(TRACE) << "SolverManager::storePreconfiguredSolvers";
  if (specificSettings_.pyObject())
  {
    if (specificSettings_.hasKey("Solvers"))
    {

      std::string keyString("Solvers");
      std::pair<std::string, PyObject *> dictItem
        = specificSettings_.getOptionDictBegin<std::string, PyObject *>(keyString);

      for (; !specificSettings_.getOptionDictEnd(keyString);
          specificSettings_.getOptionDictNext<std::string, PyObject *>(keyString, dictItem))
      {
        std::string key = dictItem.first;
        PyObject *value = dictItem.second;

        if (value == NULL)
        {
          LOG(WARNING) << "Could not extract dict for solver \"" << key << "\".";
        }
        else if(!PyDict_Check(value))
        {
          LOG(WARNING) << "Value for solver with name \"" << key << "\" should be a dict.";
        }
        else
        {
          LOG(DEBUG) << "store solver configuration with key \"" << key << "\".";
          if (solverConfiguration_.find(key) != solverConfiguration_.end())
          {
            solverConfiguration_.at(key).setPyObject(value);
          }
          else
          {
            solverConfiguration_.insert(std::pair<std::string,PythonConfig>(key, PythonConfig(specificSettings_, "Solvers", key, value)));
          }
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
