#include "field_variable/element_to_dof_mapping.h"

#include <cassert>

#include "utility/string_utility.h"
#include "utility/math_utility.h"
#include "control/types.h"
#include "easylogging++.h"

namespace FieldVariable
{
  
void ElementToDofMapping::setNumberElements(element_idx_t nElements)
{
  dofs_.resize(nElements);
}

std::shared_ptr<NodeToDofMapping> ElementToDofMapping::setup(std::shared_ptr<ExfileRepresentation> exfileRepresentation,
                                                             ElementToNodeMapping &elementToNodeMapping,
                                                             const int nDofsPerNode
                                                             )
{
  node_idx_t dofGlobalNo = 0;
  
  // for setup to work we need the number of elements already set (by a previous call to setNumberElements)
  assert(dofs_.size() != 0);
  int nElements = dofs_.size();
  
  // create node to dof mapping 
  std::shared_ptr<NodeToDofMapping> nodeToDofMapping = std::make_shared<NodeToDofMapping>();
  
  // loop over elements 
  for (element_idx_t elementNo = 0; elementNo < nElements; elementNo++)
  {
    // get available information about element
    ElementToNodeMapping::Element &element = elementToNodeMapping.getElement(elementNo);
    std::shared_ptr<ExfileElementRepresentation> exfileElement = exfileRepresentation->getExfileElementRepresentation(elementNo);
    
    // element contains the following fields:
    //   std::vector<int> globalNodeNo;
    //   std::vector<double> scaleFactors;
    
    // resize dofs vector for element
    int nNodesInElement = element.globalNodeNo.size();
    dofs_[elementNo].resize(nNodesInElement*nDofsPerNode);
    
    int elementDofIndex = 0;
    
    // loop over nodes of element
    for (int nodeIndex = 0; nodeIndex < nNodesInElement; nodeIndex++)
    {
      int nodeGlobalNo = element.globalNodeNo[nodeIndex];
      
      ExfileElementRepresentation::Node exfileNode = exfileElement->getNode(nodeIndex);
      
      // if node was not visited yet
      if (!nodeToDofMapping->containsNode(nodeGlobalNo))
      {
        NodeToDofMapping::NodeDofInformation &nodeDofInformation = nodeToDofMapping->getNodeDofInformation(nodeGlobalNo);
      
        nodeDofInformation.dofs.resize(nDofsPerNode);
        nodeDofInformation.exfileValueIndices.resize(nDofsPerNode);
        std::copy(nodeDofInformation.exfileValueIndices.begin(), exfileNode.valueIndices.begin(), exfileNode.valueIndices.end());
        
        for (int dofIndex = 0; dofIndex < nDofsPerNode; dofIndex++)
        {
          dofs_[elementNo][elementDofIndex++] = dofGlobalNo;
          nodeDofInformation.dofs[dofIndex] = dofGlobalNo++;
        }
      }
      else
      {
        // if node was already visited, i.e. it is adjacent to an element with an earlier element no.
        NodeToDofMapping::NodeDofInformation &nodeDofInformation = nodeToDofMapping->getNodeDofInformation(nodeGlobalNo);
      
        size_t subsequenceStartPos = 0;
        if (!MathUtility::isSubsequenceOf(nodeDofInformation.exfileValueIndices, exfileNode.valueIndices, subsequenceStartPos))
        {
          // node was visited earlier, but the assigned dofs are different for the current element (other version at this node)
          
          for (int dofIndex = 0; dofIndex < nDofsPerNode; dofIndex++)
          {
            dofs_[elementNo][elementDofIndex++] = dofGlobalNo;
            nodeDofInformation.dofs.push_back(dofGlobalNo++);
            nodeDofInformation.exfileValueIndices.push_back(exfileNode.valueIndices[dofIndex]);
          }
        }
        else
        {
          // node was visited earlier and has assigned dofs
          for (int dofIndex = 0; dofIndex < nDofsPerNode; dofIndex++)
          {
            int dofNo = nodeDofInformation.dofs[subsequenceStartPos+dofIndex];
            dofs_[elementNo][elementDofIndex++] = dofNo;
          }
        }
      }
    }
  }
  nDofs_ = dofGlobalNo;
  
  return nodeToDofMapping;
}

int ElementToDofMapping::nDofs() const
{
  return nDofs_;
}

std::vector<int> &ElementToDofMapping::getElementDofs(element_idx_t elementGlobalNo)
{
  return dofs_[elementGlobalNo]; 
}

};