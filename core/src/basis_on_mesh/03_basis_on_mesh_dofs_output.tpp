#include "basis_on_mesh/03_basis_on_mesh_dofs.h"

#include "easylogging++.h"
#include "utility/string_utility.h"

#include <iostream>
#include <fstream>

namespace BasisOnMesh
{
 
template<int D,typename BasisFunctionType>
void BasisOnMeshDofs<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::
outputExelemFile(std::ostream &file);
{ 
  file << " Group name: main_group" << std::endl
    << " Shape. Dimension=3, line*line*line" << std::endl;
   
  bool outputHeader = true;
   
  for(node_idx_t currentElementGlobalNo = 0; currentElementGlobalNo < nElements(); currentElementGlobalNo++)
  {
    int nScaleFactors = fieldVariable_.begin()->second->getNumberScaleFactors(currentElementGlobalNo);
    
    if (nScaleFactors != 0)
    {
      file << " #Scale factor sets=1" << std::endl
        << " c.Hermite*c.Hermite*c.Hermite, #Scale factors=" << nScaleFactors << std::endl;
    }
    
    file << " #Nodes=" << nNodesPerElement() << std::endl
      << " #Fields=" << fieldVariable_.size() << std::endl;
      
    // check if a new header is necessary
    if (currentElementGlobalNo > 0)
    {
      outputHeader = false;
      for (auto &fieldVariable : fieldVariable_)
      {
        if (!fieldVariable->second->haveSameExfileRepresentation(currentElementGlobalNo-1, currentElementGlobalNo))
        {
          outputHeader = true;
          break;
        }
      }
    }
    
    // output header 
    if (outputHeader)
    {
      for (auto &fieldVariable : fieldVariable_)
      {
        fieldVariable->second->outputHeaderExelemFile(file, currentElementGlobalNo);
      }
    }
    
    elementToNodeMapping_->outputElementExelemFile(file, currentElementGlobalNo);
  }
}

template<int D,typename BasisFunctionType>
void BasisOnMeshDofs<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::
outputExnodeFile(std::ostream &file);
{ 
  file << " Group name: main_group" << std::endl;
   
  bool outputHeader = true;
   
  for(node_idx_t currentNodeGlobalNo = 0; currentNodeGlobalNo < nElements(); currentNodeGlobalNo++)
  {
    // check if a new header is necessary
    if (currentNodeGlobalNo > 0)
    {
      outputHeader = false;
      for (auto &fieldVariable : fieldVariable_)
      {
        if (!fieldVariable->second->haveSameExfileRepresentation(currentNodeGlobalNo-1, currentNodeGlobalNo))
        {
          outputHeader = true;
          break;
        }
      }
    }
    
    // output header 
    if (outputHeader)
    {
      file << " #Fields=" << fieldVariable_.size() << std::endl;
      int valueIndex = 0;  // an index that runs over values of components of field variables and corresponds to the index in the values block for a node
      for (auto &fieldVariable : fieldVariable_)
      {
        fieldVariable->second->outputHeaderExnodeFile(file, currentNodeGlobalNo, valueIndex);
      }
    }
    
    file << " Node: " << currentNodeGlobalNo << std::endl;
    
    // collect values of all field variables at the current node
    std::vector<double> valuesAtNode;
    for (auto &fieldVariable : fieldVariable_)
    {
      fieldVariable->second->collectValuesOfNode(currentNodeGlobalNo, valuesAtNode);
    }
    StringUtility::outputValuesBlock(file, valuesAtNode, 8);
  }
}

};  // namespace