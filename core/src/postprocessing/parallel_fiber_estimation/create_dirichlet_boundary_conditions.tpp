#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
createDirichletBoundaryConditions(const std::array<int,3> &nElementsPerCoordinateDirectionLocal, std::shared_ptr<SpatialDiscretization::DirichletBoundaryConditions<FunctionSpaceType,1>> &dirichletBoundaryConditions)
{
  // create dirichlet BC object
  dirichletBoundaryConditions = std::make_shared<SpatialDiscretization::DirichletBoundaryConditions<FunctionSpaceType,1>>();

  typedef typename SpatialDiscretization::DirichletBoundaryConditions<FunctionSpaceType,1>::ElementWithNodes ElementWithNodes;
  const int nDofsPerElement1D = FunctionSpace::FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerElement();

  std::vector<ElementWithNodes> boundaryConditionElements;
  std::vector<dof_no_t> boundaryConditionNonGhostDofLocalNos;
  std::vector<std::array<double,1>> boundaryConditionValues;

  // fill dirichlet boundary condition object
  // set bottom nodes to 0
  // only if the current domain is at the bottom of the overall domain
  if (meshPartition_->ownRankPartitioningIndex(2) == 0)
  {
    std::set<dof_no_t> boundaryConditionNonGhostDofLocalNosSet;

    // loop over bottom elements
    for (element_no_t elementIndexX = 0; elementIndexX < nElementsPerCoordinateDirectionLocal[0]; elementIndexX++)
    {
      for (element_no_t elementIndexY = 0; elementIndexY < nElementsPerCoordinateDirectionLocal[1]; elementIndexY++)
      {
        element_no_t elementNoLocal = elementIndexY*nElementsPerCoordinateDirectionLocal[1] + elementIndexX;

        ElementWithNodes elementWithNodes;
        elementWithNodes.elementNoLocal = elementNoLocal;

        // loop over dofs of element that are at bottom
        for (dof_no_t elementalDofIndexX = 0; elementalDofIndexX < nDofsPerElement1D; elementalDofIndexX++)
        {
          for (dof_no_t elementalDofIndexY = 0; elementalDofIndexY < nDofsPerElement1D; elementalDofIndexY++)
          {
            dof_no_t elementalDofIndex = elementalDofIndexY*nDofsPerElement1D + elementalDofIndexX;
            elementWithNodes.elementalDofIndex.push_back(std::pair<int,std::array<double,1>>(elementalDofIndex, std::array<double,1>({0.0})));

            dof_no_t dofLocalNo = this->functionSpace_->getDofNo(elementNoLocal, elementalDofIndex);
            boundaryConditionNonGhostDofLocalNosSet.insert(dofLocalNo);
          }
        }
        boundaryConditionElements.push_back(elementWithNodes);
      }
    }

    // create the vectors for dofs and values

    // transfer the data in the set to the vector for the BC indices
    boundaryConditionNonGhostDofLocalNos.resize(boundaryConditionNonGhostDofLocalNosSet.size());
    std::copy(boundaryConditionNonGhostDofLocalNosSet.begin(), boundaryConditionNonGhostDofLocalNosSet.end(), boundaryConditionNonGhostDofLocalNos.begin());

    // set the same amount of values 0.0 for the BC values
    boundaryConditionValues.resize(boundaryConditionNonGhostDofLocalNosSet.size());
    std::fill(boundaryConditionValues.begin(), boundaryConditionValues.end(), std::array<double,1>({0.0}));
    boundaryConditionNonGhostDofLocalNosSet.clear();
  }

  // set top nodes to 1
  // loop over top elements
  // only if the current domain is at the top of the overall domain
  if (meshPartition_->ownRankPartitioningIndex(2) == meshPartition_->nRanks(2)-1)
  {
    std::set<dof_no_t> boundaryConditionNonGhostDofLocalNosSet;

    element_no_t elementIndexZ = nElementsPerCoordinateDirectionLocal[2]-1;
    for (element_no_t elementIndexX = 0; elementIndexX < nElementsPerCoordinateDirectionLocal[0]; elementIndexX++)
    {
      for (element_no_t elementIndexY = 0; elementIndexY < nElementsPerCoordinateDirectionLocal[1]; elementIndexY++)
      {
        element_no_t elementNoLocal = elementIndexZ*nElementsPerCoordinateDirectionLocal[0]*nElementsPerCoordinateDirectionLocal[1]
          + elementIndexY*nElementsPerCoordinateDirectionLocal[1] + elementIndexX;

        ElementWithNodes elementWithNodes;
        elementWithNodes.elementNoLocal = elementNoLocal;

        // loop over dofs of element that are at top
        dof_no_t elementalDofIndexZ = nDofsPerElement1D-1;
        for (dof_no_t elementalDofIndexX = 0; elementalDofIndexX < nDofsPerElement1D; elementalDofIndexX++)
        {
          for (dof_no_t elementalDofIndexY = 0; elementalDofIndexY < nDofsPerElement1D; elementalDofIndexY++)
          {
            dof_no_t elementalDofIndex = elementalDofIndexZ*nDofsPerElement1D*nDofsPerElement1D
              + elementalDofIndexY*nDofsPerElement1D + elementalDofIndexX;
            elementWithNodes.elementalDofIndex.push_back(std::pair<int,std::array<double,1>>(elementalDofIndex, std::array<double,1>({1.0})));

            dof_no_t dofLocalNo = this->functionSpace_->getDofNo(elementNoLocal, elementalDofIndex);
            boundaryConditionNonGhostDofLocalNosSet.insert(dofLocalNo);
          }
        }
        boundaryConditionElements.push_back(elementWithNodes);
      }
    }

    // fill the vectors for dofs and values
    // transfer the data in the set to the vector for the BC indices
    int nBottomDofs = boundaryConditionNonGhostDofLocalNos.size();
    boundaryConditionNonGhostDofLocalNos.resize(nBottomDofs + boundaryConditionNonGhostDofLocalNosSet.size());
    std::copy(boundaryConditionNonGhostDofLocalNosSet.begin(), boundaryConditionNonGhostDofLocalNosSet.end(), boundaryConditionNonGhostDofLocalNos.begin()+nBottomDofs);

    // add the same amount of 1.0 values for the BC values
    boundaryConditionValues.resize(nBottomDofs + boundaryConditionNonGhostDofLocalNosSet.size());
    std::fill(boundaryConditionValues.begin() + nBottomDofs, boundaryConditionValues.end(), std::array<double,1>({1.0}));
  }

  /*
  LOG(DEBUG) << "boundaryConditionElements:";
  for (auto &boundaryConditionElement : boundaryConditionElements)
  {
    LOG(DEBUG) << "{el." << boundaryConditionElement.elementNoLocal << ", (dof,v):" << boundaryConditionElement.elementalDofIndex << "}";
  }
  LOG(DEBUG) << "boundaryConditionNonGhostDofLocalNos: " << boundaryConditionNonGhostDofLocalNos;
  LOG(DEBUG) << "boundaryConditionValues: " << boundaryConditionValues;
*/
  dirichletBoundaryConditions->initialize(this->functionSpace_, boundaryConditionElements, boundaryConditionNonGhostDofLocalNos, boundaryConditionValues);
}

} // namespace
