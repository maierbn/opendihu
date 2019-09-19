#include "output_writer/exfile/exfile_writer.h"

#include "function_space/function_space.h"

#include <cstdlib>

namespace OutputWriter
{

//! write exelem file to given stream
template<int D, typename BasisFunctionType, typename FieldVariablesForOutputWriterType>
void ExfileWriter<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,FieldVariablesForOutputWriterType>::
outputExelem(std::ostream &stream, FieldVariablesForOutputWriterType fieldVariables, std::string meshName, 
             std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>> mesh, 
             int nFieldVariablesOfMesh)
{
  stream << " Group name: " << meshName << std::endl
    << " Shape. Dimension=" << D << ", " << StringUtility::multiply<D>("line") << std::endl;

  const int nNodesPerElement = FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::nNodesPerElement();
  const element_no_t nElements = mesh->nElementsLocal();

  bool outputHeader = true;

  // loop over elements
  for (element_no_t currentElementGlobalNo = 0; currentElementGlobalNo < nElements; currentElementGlobalNo++)
  {
    // check if a new header is necessary
    if (currentElementGlobalNo > 0)
    {
      outputHeader = false;
      ExfileLoopOverTuple::loopCheckIfNewExelemHeaderNecessary<FieldVariablesForOutputWriterType>(fieldVariables, meshName, currentElementGlobalNo, outputHeader);
    }

    // output header
    if (outputHeader)
    {
      int nScaleFactors = mesh->getNumberScaleFactors(currentElementGlobalNo);
      
      if (nScaleFactors == 0)
      {
        stream << " #Scale factor sets=0" << std::endl;
      }
      else
      {
        stream << " #Scale factor sets=1" << std::endl;
        stream << " " << BasisFunction::getBasisRepresentationString<D,BasisFunctionType>();
        stream << ", #Scale factors=" << nScaleFactors << std::endl;
      }

      stream << " #Nodes=" << nNodesPerElement << std::endl
        << " #Fields=" << nFieldVariablesOfMesh << std::endl;

      // loop over field variables and output headers
      int fieldVariableIndex = 0;
      ExfileLoopOverTuple::loopOutputHeaderExelem<FieldVariablesForOutputWriterType>(fieldVariables, fieldVariableIndex, meshName, stream, currentElementGlobalNo);
    }
    
    mesh->elementToNodeMapping()->outputElementExelem(stream, currentElementGlobalNo);
  }

}

//! write exnode file to given stream
template<int D, typename BasisFunctionType, typename FieldVariablesForOutputWriterType>
void ExfileWriter<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>,FieldVariablesForOutputWriterType>::
outputExnode(std::ostream &stream, FieldVariablesForOutputWriterType fieldVariables, std::string meshName, 
             std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>> mesh,
             int nFieldVariablesOfMesh
            )
{
  stream << " Group name: " << meshName << std::endl;

  bool outputHeader = true;
  const node_no_t nNodes = mesh->nNodesGlobal();

  // loop over all nodes
  for (node_no_t currentNodeGlobalNo = 0; currentNodeGlobalNo < nNodes; currentNodeGlobalNo++)
  {
    // check if a new header is necessary
    if (currentNodeGlobalNo > 0)
    {
      outputHeader = false;
      ExfileLoopOverTuple::loopCheckIfNewExnodeHeaderNecessary<FieldVariablesForOutputWriterType>(fieldVariables, meshName, currentNodeGlobalNo, outputHeader);
    }

    // output header
    if (outputHeader)
    {
      stream << " #Fields=" << nFieldVariablesOfMesh << std::endl;
      int valueIndex = 0;  // an index that runs over values of components of field variables and corresponds to the index in the values block for a node

      // output the exnode file header for the current node
      int fieldVariableIndex = 0;
      ExfileLoopOverTuple::loopOutputHeaderExnode<FieldVariablesForOutputWriterType>(fieldVariables, fieldVariableIndex, meshName, stream, currentNodeGlobalNo, valueIndex);
    }

    stream << " Node: " << currentNodeGlobalNo+1 << std::endl;

    // collect values of all field variables at the current node
    // get dofs
    std::vector<double> valuesAtNode;
    ExfileLoopOverTuple::loopGetValuesAtNode<FieldVariablesForOutputWriterType>(fieldVariables, meshName, currentNodeGlobalNo, valuesAtNode);

    StringUtility::outputValuesBlock(stream, valuesAtNode.begin(), valuesAtNode.end(), 8);
  }
}


}  // namespace
