#include "output_writer/exfile_writer.h"

#include "basis_on_mesh/05_basis_on_mesh.h"

namespace OutputWriter
{

//! write exelem file to given stream
template<int D, typename BasisFunctionType>
void ExfileWriter<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
outputExelem(std::ostream &stream, std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> fieldVariables)
{
  stream << " Group name: main_group" << std::endl
    << " Shape. Dimension=" << D << ", " << StringUtility::multiply<D>("line") << std::endl;
  
  const int nNodesPerElement = BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>::nNodesPerElement();   
  const element_no_t nElements = fieldVariables.front()->mesh()->nElements();
 
  bool outputHeader = true;
   
  for(element_no_t currentElementGlobalNo = 0; currentElementGlobalNo < nElements; currentElementGlobalNo++)
  {
    // check if a new header is necessary
    if (currentElementGlobalNo > 0)
    {
      outputHeader = false;
      for (auto &fieldVariable : fieldVariables)
      {
        VLOG(2) << "check if field variable " << fieldVariable->name() << " has the same exfileRepr for elements " 
          << currentElementGlobalNo-1 << " and " << currentElementGlobalNo;
          
        if (!fieldVariable->haveSameExfileRepresentation(currentElementGlobalNo-1, currentElementGlobalNo))
        {
          outputHeader = true;
          break;
        }
      }
    }
    
    // output header 
    if (outputHeader)
    {    
      int nScaleFactors = fieldVariables.front()->getNumberScaleFactors(currentElementGlobalNo);
      
      if (nScaleFactors != 0)
      {
        stream << " #Scale factor sets=1" << std::endl;
        stream << " " << BasisFunction::getBasisRepresentationString<D,BasisFunctionType>();
        stream << ", #Scale factors=" << nScaleFactors << std::endl;
      }
      
      stream << " #Nodes=" << nNodesPerElement << std::endl
        << " #Fields=" << fieldVariables.size() << std::endl;
      
      int fieldVariableNo = 0;     // a number that runs over the field variables
      for (auto &fieldVariable : fieldVariables)
      {
        fieldVariable->outputHeaderExelem(stream, currentElementGlobalNo, fieldVariableNo++);
      }
    }
    
    fieldVariables.front()->elementToNodeMapping()->outputElementExelemFile(stream, currentElementGlobalNo);
  }
  
}

//! write exnode file to given stream
template<int D, typename BasisFunctionType>
void ExfileWriter<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::
outputExnode(std::ostream &stream, std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> fieldVariables)
{
  stream << " Group name: main_group" << std::endl;

  bool outputHeader = true;
   
  const node_no_t nNodes = fieldVariables.front()->mesh()->nNodes();
  
  // loop over all nodes
  for(node_no_t currentNodeGlobalNo = 0; currentNodeGlobalNo < nNodes; currentNodeGlobalNo++)
  {
    // check if a new header is necessary
    if (currentNodeGlobalNo > 0)
    {
      outputHeader = false;
      for (auto &fieldVariable : fieldVariables)
      {
        int previousNumberVersions = fieldVariable->nodeToDofMapping()->nVersions(currentNodeGlobalNo-1);
        int currentNumberVersions = fieldVariable->nodeToDofMapping()->nVersions(currentNodeGlobalNo);
        
        if(previousNumberVersions != currentNumberVersions)
        {
          outputHeader = true;
          break;
        }
      }
    }
    
    // output header 
    if (outputHeader)
    {
      stream << " #Fields=" << fieldVariables.size() << std::endl;
      int fieldVariableNo = 0;     // a number that runs over the field variables
      int valueIndex = 0;  // an index that runs over values of components of field variables and corresponds to the index in the values block for a node
      for (auto &fieldVariable : fieldVariables)
      {
        fieldVariable->outputHeaderExnode(stream, currentNodeGlobalNo, valueIndex, fieldVariableNo++);
      }
    }
    
    stream << " Node: " << currentNodeGlobalNo+1 << std::endl;
    
    // collect values of all field variables at the current node
    // get dofs
    std::vector<double> valuesAtNode;
    int i=0;
    for (auto &fieldVariable : fieldVariables)
    {
      VLOG(2) << "  fieldVariable " << i++ << ": " << fieldVariable->name() << " nodeToDofMapping: " << fieldVariable->nodeToDofMapping();
      
      std::vector<int> &dofsAtNode = fieldVariable->nodeToDofMapping()->getNodeDofs(currentNodeGlobalNo);
      
      // loop over components
      for (std::string component : fieldVariable->componentNames())
      {
        
        //stream << std::endl << "fieldVariable " << fieldVariable->name() << "." << component<<", dofsAtNode: " << dofsAtNode << std::endl;
      
        // loop over dofs 
        for (dof_no_t dofGlobalNo : dofsAtNode)
        {
          double value = fieldVariable->getValue(component, dofGlobalNo);
          //stream << std::endl << "dof " << dofGlobalNo << ", value: " << value << std::endl;
          valuesAtNode.push_back(value);
        }
      }
    }
    
    StringUtility::outputValuesBlock(stream, valuesAtNode, 8);
  }
}

  
};  //namespace
