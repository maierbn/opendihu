#include "output_writer/exfile/exfile_writer.h"

#include "basis_on_mesh/05_basis_on_mesh.h"

#include <cstdlib>

namespace OutputWriter
{

namespace ExfileLoopOverTuple
{

template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i < std::tuple_size<OutputFieldVariablesType>::value, void>::type
checkIfNewExelemHeaderNecessary(const OutputFieldVariablesType &fieldVariables, element_no_t currentElementGlobalNo, bool &newHeaderNecessary)
{
  VLOG(2) << "check if field variable " << std::get<i>(fieldVariables)->name() << " has the same exfileRepr for elements " 
    << currentElementGlobalNo-1 << " and " << currentElementGlobalNo;
   
  if (!std::get<i>(fieldVariables)->haveSameExfileRepresentation(currentElementGlobalNo-1, currentElementGlobalNo))
  {
    newHeaderNecessary = true;
    return;
  }
    
  // advance to next tuple element
  checkIfNewExelemHeaderNecessary<OutputFieldVariablesType, i+1>(fieldVariables, currentElementGlobalNo, newHeaderNecessary);
}

template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i < std::tuple_size<OutputFieldVariablesType>::value, void>::type
checkIfNewExnodeHeaderNecessary(const OutputFieldVariablesType &fieldVariables, element_no_t currentNodeGlobalNo, bool &newHeaderNecessary)
{
  int previousNumberVersions = std::get<i>(fieldVariables)->nodeToDofMapping()->nVersions(currentNodeGlobalNo-1);
  int currentNumberVersions = std::get<i>(fieldVariables)->nodeToDofMapping()->nVersions(currentNodeGlobalNo);
  
  if(previousNumberVersions != currentNumberVersions)
  {
    newHeaderNecessary = true;
    return;
  }
    
  // advance to next tuple element
  checkIfNewExnodeHeaderNecessary<OutputFieldVariablesType, i+1>(fieldVariables, currentNodeGlobalNo, newHeaderNecessary);
}

template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i < std::tuple_size<OutputFieldVariablesType>::value, void>::type
getValuesAtNode(const OutputFieldVariablesType &fieldVariables, element_no_t currentNodeGlobalNo, std::vector<double> &valuesAtNode)
{
  
  VLOG(1) << "  fieldVariable " << i << ": " << std::get<i>(fieldVariables)->name() 
    << " nodeToDofMapping: " << std::get<i>(fieldVariables)->nodeToDofMapping();
    
  std::vector<dof_no_t> &dofsAtNode = std::get<i>(fieldVariables)->nodeToDofMapping()->getNodeDofs(currentNodeGlobalNo);
    
  // loop over components
  for (int componentNo = 0; componentNo < std::get<i>(fieldVariables)->getNComponents(); componentNo++)
  {
    std::get<i>(fieldVariables)->getValues(componentNo, dofsAtNode, valuesAtNode);
    VLOG(1) << "   Component " << componentNo << ", dofsAtNode: " << dofsAtNode << ", valuesAtNode: " << valuesAtNode;
  }
   
  // advance to next tuple element
  getValuesAtNode<OutputFieldVariablesType, i+1>(fieldVariables, currentNodeGlobalNo, valuesAtNode);
}

};  // namespace

//! write exelem file to given stream
template<int D, typename BasisFunctionType, typename OutputFieldVariablesType>
void ExfileWriter<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,OutputFieldVariablesType>::
outputExelem(std::ostream &stream, OutputFieldVariablesType fieldVariables)
{
  stream << " Group name: Region" << std::endl
    << " Shape. Dimension=" << D << ", " << StringUtility::multiply<D>("line") << std::endl;
  
  const int nNodesPerElement = BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::nNodesPerElement();   
  const element_no_t nElements = std::get<0>(fieldVariables)->mesh()->nElements();
 
  bool outputHeader = true;
   
  for(element_no_t currentElementGlobalNo = 0; currentElementGlobalNo < nElements; currentElementGlobalNo++)
  {
    // check if a new header is necessary
    if (currentElementGlobalNo > 0)
    {
      outputHeader = false;
      ExfileLoopOverTuple::checkIfNewExelemHeaderNecessary<OutputFieldVariablesType>(fieldVariables, currentElementGlobalNo, outputHeader);
    }
    
    // output header 
    if (outputHeader)
    {    
      int nScaleFactors = std::get<0>(fieldVariables)->getNumberScaleFactors(currentElementGlobalNo);
      
      if (nScaleFactors != 0)
      {
        stream << " #Scale factor sets=1" << std::endl;
        stream << " " << BasisFunction::getBasisRepresentationString<D,BasisFunctionType>();
        stream << ", #Scale factors=" << nScaleFactors << std::endl;
      }
      
      stream << " #Nodes=" << nNodesPerElement << std::endl
        << " #Fields=" << std::tuple_size<OutputFieldVariablesType>::value << std::endl;
      
      // loop over field variables and output headers
      ExfileLoopOverTuple::outputHeaderExelem<OutputFieldVariablesType>(fieldVariables, stream, currentElementGlobalNo);
    }
    
    std::get<0>(fieldVariables)->elementToNodeMapping()->outputElementExelem(stream, currentElementGlobalNo);
  }
  
}

//! write exnode file to given stream
template<int D, typename BasisFunctionType, typename OutputFieldVariablesType>
void ExfileWriter<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,OutputFieldVariablesType>::
outputExnode(std::ostream &stream, OutputFieldVariablesType fieldVariables)
{
  stream << " Group name: Region" << std::endl;

  bool outputHeader = true;
  const node_no_t nNodes = std::get<0>(fieldVariables)->mesh()->nNodes();
  
  // loop over all nodes
  for(node_no_t currentNodeGlobalNo = 0; currentNodeGlobalNo < nNodes; currentNodeGlobalNo++)
  {
    // check if a new header is necessary
    if (currentNodeGlobalNo > 0)
    {
      outputHeader = false;
      ExfileLoopOverTuple::checkIfNewExnodeHeaderNecessary<OutputFieldVariablesType>(fieldVariables, currentNodeGlobalNo, outputHeader);
    }
    
    // output header 
    if (outputHeader)
    {
      stream << " #Fields=" << std::tuple_size<OutputFieldVariablesType>::value << std::endl;
      int valueIndex = 0;  // an index that runs over values of components of field variables and corresponds to the index in the values block for a node
      
      // output the exnode file header for the current node
      ExfileLoopOverTuple::outputHeaderExnode<OutputFieldVariablesType>(fieldVariables, stream, currentNodeGlobalNo, valueIndex);
    }
    
    stream << " Node: " << currentNodeGlobalNo+1 << std::endl;
    
    // collect values of all field variables at the current node
    // get dofs
    std::vector<double> valuesAtNode;
    ExfileLoopOverTuple::getValuesAtNode<OutputFieldVariablesType>(fieldVariables, currentNodeGlobalNo, valuesAtNode);
    
    StringUtility::outputValuesBlock(stream, valuesAtNode.begin(), valuesAtNode.end(), 8);
  }
}

  
};  //namespace
