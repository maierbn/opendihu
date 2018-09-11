#include "utility/boundary_conditions.h"

#include "easylogging++.h"
#include "utility/python_utility.h"
#include "control/types.h"

namespace BoundaryConditions
{

template<typename FunctionSpaceType, typename T>
void parseBoundaryConditions(PyObject *settings, std::shared_ptr<FunctionSpaceType> functionSpace,
                             std::vector<std::pair<int,T>> &boundaryConditions)
{
  LOG(TRACE) << "parseBoundaryConditions";
  assert(functionSpace);

  // add weak form of Dirichlet BC to rhs
  const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();

  // determine if the BC indices in the config are given for global or local dof nos
  bool inputMeshIsGlobal = PythonUtility::getOptionBool(settings, "inputMeshIsGlobal", true);

  int nDofs = 0;
  if (inputMeshIsGlobal)
  {
    nDofs = functionSpace->nDofsGlobal();
  }
  else
  {
    nDofs = functionSpace->nDofsLocalWithoutGhosts();
  }

  if (PythonUtility::hasKey(settings, "DirichletBoundaryCondition"))
  {
    LOG(ERROR) << "Option \"DirichletBoundaryCondition\" was renamed to \"dirichletBoundaryConditions\".";
  }

  // parse all boundary conditions that are given in config
  std::pair<int, T> boundaryCondition = PythonUtility::getOptionDictBegin<int, T>(settings, "dirichletBoundaryConditions");
  for (; !PythonUtility::getOptionDictEnd(settings, "dirichletBoundaryConditions");
          PythonUtility::getOptionDictNext<int, T>(settings, "dirichletBoundaryConditions", boundaryCondition))
  {
    // for negative indices add number of dofs such that -1 is the last dof, -2 is the econd-last etc.
    if (boundaryCondition.first < 0)
    {
      boundaryCondition.first += nDofs;
    }
    else if (boundaryCondition.first > functionSpace->nDofsLocalWithoutGhosts())
    {
      node_no_t nodeNoLocal = boundaryCondition.first / nDofsPerNode;
      node_no_t nNodesLocal = functionSpace->nNodesLocalWithoutGhosts();
      LOG(ERROR) << "Boundary condition specified for index " << boundaryCondition.first << " (node " << nodeNoLocal << "), "
        << "but there are only " << functionSpace->nDofsLocalWithoutGhosts() << " local dofs (" << nNodesLocal << " nodes)";
    }
    boundaryConditions.push_back(boundaryCondition);
  }
  if (boundaryConditions.empty())
  {
    LOG(DEBUG) << "no boundary conditions specified";
  }
}

}  // namespace
