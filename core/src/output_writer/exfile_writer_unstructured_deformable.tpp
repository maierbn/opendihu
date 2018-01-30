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
  const int nElements = fieldVariables.front()->mesh()->nElements();
 
  bool outputHeader = true;
   
  for(element_idx_t currentElementGlobalNo = 0; currentElementGlobalNo < nElements; currentElementGlobalNo++)
  {
    // check if a new header is necessary
    if (currentElementGlobalNo > 0)
    {
      outputHeader = false;
      for (auto &fieldVariable : fieldVariables)
      {
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
        if (typeid(BasisFunctionType) == typeid(typename BasisFunction::Hermite))
          stream << " " << StringUtility::multiply<D>("c.Hermite");
        else if (typeid(BasisFunctionType) == typeid(typename BasisFunction::Lagrange<1>))
          stream << " " << StringUtility::multiply<D>("l.Lagrange");
        else if (typeid(BasisFunctionType) == typeid(typename BasisFunction::Lagrange<2>))
          stream << " " << StringUtility::multiply<D>("q.Lagrange");
        else 
          LOG(FATAL) << "Exfile output not implemented for basis function type";
        
        stream << ", #Scale factors=" << nScaleFactors << std::endl;
      }
      
      stream << " #Nodes=" << nNodesPerElement << std::endl
        << " #Fields=" << fieldVariables.size() << std::endl;
      
      int fieldVariableNo = 0;
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
   
  const int nNodes = fieldVariables.front()->mesh()->nNodes();
  
  // loop over all nodes
  for(node_idx_t currentNodeGlobalNo = 0; currentNodeGlobalNo < nNodes; currentNodeGlobalNo++)
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
      int valueIndex = 0;  // an index that runs over values of components of field variables and corresponds to the index in the values block for a node
      for (auto &fieldVariable : fieldVariables)
      {
        fieldVariable->outputHeaderExnode(stream, currentNodeGlobalNo, valueIndex);
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
        for (dof_idx_t dofGlobalNo : dofsAtNode)
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
