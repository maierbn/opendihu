#include "spatial_discretization/boundary_conditions/boundary_conditions_base.h"

#include "easylogging++.h"
#include "utility/python_utility.h"
#include "utility/vector_operators.h"
#include "control/types.h"

namespace SpatialDiscretization
{


template<typename FunctionSpaceType,int nComponents>
BoundaryConditionsBase<FunctionSpaceType,nComponents>::
BoundaryConditionsBase(DihuContext context) : specificSettings_(NULL)
{
}

template<typename FunctionSpaceType,int nComponents>
void BoundaryConditionsBase<FunctionSpaceType,nComponents>::
initialize(PythonConfig specificSettings, std::shared_ptr<FunctionSpaceType> functionSpace, std::string boundaryConditionsConfigKey)
{
  functionSpace_ = functionSpace;
  specificSettings_ = specificSettings;
  printDebuggingInfo();

  parseBoundaryConditionsForElements(boundaryConditionsConfigKey);
  this->initializeGhostElements();
}

template<typename FunctionSpaceType,int nComponents>
void BoundaryConditionsBase<FunctionSpaceType,nComponents>::
initialize(std::shared_ptr<FunctionSpaceType> functionSpace, std::vector<BoundaryConditionsBase<FunctionSpaceType,nComponents>::ElementWithNodes> &boundaryConditionElements,
                std::vector<dof_no_t> &boundaryConditionNonGhostDofLocalNos, std::vector<ValueType> &boundaryConditionValues)
{
  functionSpace_ = functionSpace;

  boundaryConditionElements_ = boundaryConditionElements;
  boundaryConditionNonGhostDofLocalNos_ = boundaryConditionNonGhostDofLocalNos;
  boundaryConditionValues_ = boundaryConditionValues;
  generateBoundaryConditionsByComponent();
  printDebuggingInfo();

  this->initializeGhostElements();
}

template<typename FunctionSpaceType,int nComponents>
void BoundaryConditionsBase<FunctionSpaceType,nComponents>::
printDebuggingInfo()
{
  if (VLOG_IS_ON(1))
  {
    VLOG(1) << "parsed boundary conditions:";
    VLOG(1) << "  dofsLocal of BC: " << boundaryConditionNonGhostDofLocalNos_ << " with prescribed values: " << boundaryConditionValues_;

    for (auto boundaryConditionElement: boundaryConditionElements_)
    {
      VLOG(1) << "  elementNo: " << boundaryConditionElement.elementNoLocal << " has (dof,value): " << boundaryConditionElement.elementalDofIndex;
    }
  }
}

template<typename FunctionSpaceType,int nComponents>
void BoundaryConditionsBase<FunctionSpaceType,nComponents>::
parseBoundaryConditions(PythonConfig settings, std::shared_ptr<FunctionSpaceType> functionSpace, std::string boundaryConditionsConfigKey,
                        std::vector<std::pair<int,std::array<double,nComponents>>> &boundaryConditions)
{
  LOG(TRACE) << "parseBoundaryConditions, key \"" << boundaryConditionsConfigKey << "\"";
  assert(functionSpace);

  if (VLOG_IS_ON(1))
  {
    PythonUtility::printDict(settings.pyObject());
  }

  const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();

  // determine if the BC indices in the config are given for global or local dof nos
  bool inputMeshIsGlobal = settings.getOptionBool("inputMeshIsGlobal", true);

  int nDofs = 0;
  if (inputMeshIsGlobal)
  {
    nDofs = functionSpace->nDofsGlobal();
  }
  else
  {
    nDofs = functionSpace->nDofsLocalWithoutGhosts();
  }

  if (settings.hasKey("DirichletBoundaryCondition"))
  {
    LOG(ERROR) << "Option \"DirichletBoundaryCondition\" was renamed to \"dirichletBoundaryConditions\".";
  }

  VLOG(1) << "inputMeshIsGlobal: " << inputMeshIsGlobal << ", nDofs: " << nDofs;

  // parse all boundary conditions that are given in config
  std::pair<int,std::array<double,nComponents>> boundaryCondition = settings.getOptionDictBegin<int,std::array<double,nComponents>>(boundaryConditionsConfigKey);
  for (; !settings.getOptionDictEnd(boundaryConditionsConfigKey);
          settings.getOptionDictNext<int,std::array<double,nComponents>>(boundaryConditionsConfigKey, boundaryCondition))
  {
    // for negative indices add number of dofs such that -1 is the last dof, -2 is the second-last etc.
    if (boundaryCondition.first < 0)
    {
      boundaryCondition.first += nDofs;
    }
    else if (boundaryCondition.first > nDofs)
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

template<typename FunctionSpaceType,int nComponents>
void BoundaryConditionsBase<FunctionSpaceType,nComponents>::
parseBoundaryConditionsForElements(std::string boundaryConditionsConfigKey)
{
  LOG(TRACE) << "parseBoundaryConditionsForElements";

  // add weak form of BC to rhs
  const int nDofsPerNode = FunctionSpaceType::nDofsPerNode();

  // determine if the BC indices in the config are given for global or local dof nos
  bool inputMeshIsGlobal = this->specificSettings_.getOptionBool("inputMeshIsGlobal", true);

  // Boundary conditions are specified for dof numbers, not nodes, such that for Hermite it is possible to prescribe derivatives.
  // However the ordering of the dofs is not known in the config for unstructured meshes. Therefore the ordering is special.
  // For every node there are as many values as dofs, in contiguous order.
  // Example for 2D Hermite, unstructured grid, 2x2 elements:
  //
  // node numbering:
  //  6_7_8
  // 3|_4_|5
  // 0|_1_|2
  //
  // dof numbering:
  //  6_7_8
  // 2|_3_|5
  // 0|_1_|4
  //
  // To specify du/dn = 0 at the left boundary you would set:
  // bc[0*2+1] = 0, bc[3*2+1] = 0, bc[6*2+1] = 0
  //
  // To specifiy u=0 on the bottom, you would set:
  // bc[0] = 0, bc[2] = 0, bc[4] = 0

  // ValueType is the array with number of components: std::array<double,nComponents>
  std::vector<std::pair<int,ValueType>> boundaryConditions;  // (index, value)
  parseBoundaryConditions(this->specificSettings_, functionSpace_, boundaryConditionsConfigKey, boundaryConditions);

  LOG(DEBUG) << "read in Dirichlet boundary conditions from config: " << boundaryConditions;

  // sort all parsed boundary conditions for their index no
  auto compareFunction = [](const std::pair<int,ValueType> &item1, const std::pair<int,ValueType> &item2)
  {
    return item1.first < item2.first;
  };
  std::sort(boundaryConditions.begin(), boundaryConditions.end(), compareFunction);

  LOG(DEBUG) << "sorted Dirichlet boundary conditions from config: " << boundaryConditions;

  // determine elements with nodes that have prescribed boundary conditions, store them in the vector boundaryConditionElements_,
  // which is organized by local elements
  element_no_t lastBoundaryConditionElement = -1;
  std::set<dof_no_t> boundaryConditionNonGhostDofLocalNosSet;   ///< same data as in boundaryConditionNonGhostDofLocalNos_, but as set

  // loop over all local elements
  for (element_no_t elementNoLocal = 0; elementNoLocal < functionSpace_->nElementsLocal(); elementNoLocal++)
  {
    VLOG(2) << "elementNoLocal: " << elementNoLocal;

    // loop over the nodes of this element
    int elementalDofIndex = 0;
    for (int nodeIndex = 0; nodeIndex < FunctionSpaceType::nNodesPerElement(); nodeIndex++)
    {
      // get global or local nodeNo of current node (depending on inputMeshIsGlobal)
      global_no_t nodeNo = 0;
      if (inputMeshIsGlobal)
      {
        global_no_t elementNoGlobalNatural = functionSpace_->meshPartition()->getElementNoGlobalNatural(elementNoLocal);
        nodeNo = functionSpace_->getNodeNoGlobalNatural(elementNoGlobalNatural, nodeIndex);
      }
      else
      {
        nodeNo = functionSpace_->getNodeNo(elementNoLocal, nodeIndex);
      }

      // loop over dofs of node and thus over the elemental dofs
      for (int nodalDofIndex = 0; nodalDofIndex < nDofsPerNode; nodalDofIndex++, elementalDofIndex++)
      {
        global_no_t indexForBoundaryCondition = nodeNo*nDofsPerNode + nodalDofIndex;

        // check if a boundary condition for the current node is given
        bool boundaryConditionForDofFound = false;

        // loop over all stored boundary conditions
        ValueType boundaryConditionValue;
        for (typename std::vector<std::pair<int,ValueType>>::const_iterator boundaryConditionIter = boundaryConditions.begin();
             boundaryConditionIter != boundaryConditions.end(); boundaryConditionIter++)
        {
          if (boundaryConditionIter->first == indexForBoundaryCondition)
          {
            boundaryConditionValue = boundaryConditionIter->second;
            boundaryConditionForDofFound = true;
            break;
          }
          if (boundaryConditionIter->first >= indexForBoundaryCondition)
            break;
        }

        VLOG(2) << "boundaryConditionForDofFound: " << boundaryConditionForDofFound;

        // here boundaryConditionIter->first is equal to indexForBoundaryCondition (then the current element/node/dof matches the boundary condition)
        // or boundaryConditionIter->first > indexForBoundaryCondition then the current dofIndex does not have any boundary condition
        // or boundaryConditionIter == boundaryConditions.end() then, too

        // if the currently considered boundaryCondition entry from config matches the current nodeNo and nodalDofIndex
        if (boundaryConditionForDofFound)
        {
          VLOG(2) << "elementNoLocal: " << elementNoLocal << ", lastBoundaryConditionElement: " << lastBoundaryConditionElement;

          // if there is not yet an entry in boundaryConditionElements with the current element, create one
          if (elementNoLocal != lastBoundaryConditionElement)
          {
            // add current element
            boundaryConditionElements_.emplace_back();
            boundaryConditionElements_.back().elementNoLocal = elementNoLocal;

            lastBoundaryConditionElement = elementNoLocal;

            VLOG(2) << "add empty entry for elementNoLocal " << elementNoLocal;
          }

          // add current node and boundary condition value to list of boundary conditions for current element
          ElementWithNodes &boundaryConditionElement = boundaryConditionElements_.back();

          VLOG(2) << " add (el-dof, value)" << std::pair<int,ValueType>(elementalDofIndex, boundaryConditionValue) << ", to boundaryConditionElement of element " << boundaryConditionElement.elementNoLocal;
            boundaryConditionElement.elementalDofIndex.push_back(std::pair<int,ValueType>(elementalDofIndex, boundaryConditionValue));

          // also store localDofNo
          dof_no_t dofLocalNo = functionSpace_->getDofNo(elementNoLocal, elementalDofIndex);

          if (dofLocalNo < functionSpace_->nDofsLocalWithoutGhosts())
          {
            // if dofLocalNo is not already contained in boundaryConditionNonGhostDofLocalNos_
            if (boundaryConditionNonGhostDofLocalNosSet.find(dofLocalNo) == boundaryConditionNonGhostDofLocalNosSet.end())
            {
              boundaryConditionNonGhostDofLocalNosSet.insert(dofLocalNo);
              boundaryConditionNonGhostDofLocalNos_.push_back(dofLocalNo);
              boundaryConditionValues_.push_back(boundaryConditionValue);
            }
          }
        }
      }
    }
  }

  generateBoundaryConditionsByComponent();

  LOG(DEBUG) << "boundaryConditionNonGhostDofLocalNos: " << boundaryConditionNonGhostDofLocalNos_
    << ", boundaryConditionValues: " << boundaryConditionValues_ << ", boundaryConditionsByComponent_: ";
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    LOG(DEBUG) << "  component " << componentNo << ", dofs " << boundaryConditionsByComponent_[componentNo].dofNosLocal
      << ", values " << boundaryConditionsByComponent_[componentNo].values;
  }
}

template<typename FunctionSpaceType,int nComponents>
void BoundaryConditionsBase<FunctionSpaceType,nComponents>::
generateBoundaryConditionsByComponent()
{
  // fill the data structure
  //  struct BoundaryConditionsForComponent
  //  {
  //    std::vector<dof_no_t> dofsLocal;
  //    std::vector<double> values;
  //  };
  //  std::array<BoundaryConditionsForComponent, nComponents> boundaryConditionsByComponent_;   ///< the boundary condition data organized by component

  // from the two data structures
  //   std::vector<dof_no_t> boundaryConditionNonGhostDofLocalNos_;        ///< vector of all local (non-ghost) boundary condition dofs
  //   std::vector<ValueType> boundaryConditionValues_;               ///< vector of the local prescribed values, related to boundaryConditionNonGhostDofLocalNos_

  // The boundaryConditionValues_ vector can contain NaN values, where in the config None was set.
  // This indicates that these dofs should not be set as Dirichlet boundary conditions.

  if (nComponents == 1)
  {
    boundaryConditionsByComponent_[0].dofNosLocal = boundaryConditionNonGhostDofLocalNos_;
    boundaryConditionsByComponent_[0].values.resize(boundaryConditionValues_.size());
    for (int i = 0; i < boundaryConditionValues_.size(); i++)
    {
      boundaryConditionsByComponent_[0].values[i] = boundaryConditionValues_[i][0];
    }
  }
  else
  {
    for (int componentNo = 0; componentNo < nComponents; componentNo++)
    {
      // preallocate as if all values were valid boundary conditions
      boundaryConditionsByComponent_[componentNo].dofNosLocal.reserve(boundaryConditionValues_.size());
      boundaryConditionsByComponent_[componentNo].values.reserve(boundaryConditionValues_.size());

      // iterate over all values
      for (int i = 0; i < boundaryConditionValues_.size(); i++)
      {
        // if value is not nan and therefore valid
        if (std::isfinite(boundaryConditionValues_[i][componentNo]))
        {
          boundaryConditionsByComponent_[componentNo].dofNosLocal.push_back(boundaryConditionNonGhostDofLocalNos_[i]);
          boundaryConditionsByComponent_[componentNo].values.push_back(boundaryConditionValues_[i][componentNo]);
        }
      }
    }
  }
}

//! get a reference to the vector of bc local dofs
template<typename FunctionSpaceType,int nComponents>
const std::vector<dof_no_t> &BoundaryConditionsBase<FunctionSpaceType,nComponents>::
boundaryConditionNonGhostDofLocalNos() const
{
  return boundaryConditionNonGhostDofLocalNos_;
}

//! get a reference to the vector of bc local dofs
template<typename FunctionSpaceType,int nComponents>
const std::vector<typename BoundaryConditionsBase<FunctionSpaceType,nComponents>::ValueType> &BoundaryConditionsBase<FunctionSpaceType,nComponents>::
boundaryConditionValues() const
{
  return boundaryConditionValues_;
}

template<typename FunctionSpaceType, int nComponents>
std::ostream &operator<<(std::ostream &stream, const typename BoundaryConditionsBase<FunctionSpaceType,nComponents>::BoundaryConditionsForComponent rhs)
{
  stream << "{dofNos: " << rhs.dofNosLocal << ", values: " << rhs.values << "}";
  return stream;
}

template<typename FunctionSpaceType, int nComponents>
std::ostream &operator<<(std::ostream &stream, const typename BoundaryConditionsBase<FunctionSpaceType,nComponents>::ElementWithNodes rhs)
{
  stream << "{el." << rhs.elementNoLocal << ", (dof,v):" << rhs.elementalDofIndex << "}";
  return stream;
}

}  // namespace
