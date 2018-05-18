#include "output_writer/exfile/exfile_writer.h"

#include <iomanip>

namespace OutputWriter
{   

namespace ExfileLoopOverTuple
{

template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i < std::tuple_size<OutputFieldVariablesType>::value, void>::type
outputHeaderExelem(const OutputFieldVariablesType &fieldVariables, std::ostream &stream, element_no_t currentElementGlobalNo)
{
  std::get<i>(fieldVariables)->outputHeaderExelem(stream, currentElementGlobalNo, i);
  
  // advance to next tuple element
  outputHeaderExelem<OutputFieldVariablesType, i+1>(fieldVariables, stream, currentElementGlobalNo);
}

template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i < std::tuple_size<OutputFieldVariablesType>::value, void>::type
outputHeaderExnode(const OutputFieldVariablesType &fieldVariables, std::ostream &stream, node_no_t currentNodeGlobalNo, int &valueIndex)
{
  std::get<i>(fieldVariables)->outputHeaderExnode(stream, currentNodeGlobalNo, valueIndex, i);
  
  // advance to next tuple element
  outputHeaderExnode<OutputFieldVariablesType, i+1>(fieldVariables, stream, currentNodeGlobalNo, valueIndex);
}

template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i < std::tuple_size<OutputFieldVariablesType>::value, void>::type
outputNodeValues(const OutputFieldVariablesType &fieldVariables, std::ostream &stream, node_no_t nodeGlobalNo)
{
  // output the values of a node for the current field variable
 
  // get the class of the current field variable
  typedef typename std::tuple_element<i,OutputFieldVariablesType>::type FieldVariableType;
  
  // get the current field variable
  auto fieldVariable = std::get<i>(fieldVariables);  // this type is std::shared_ptr<FieldVariable<..>>
  
  const int nDofsPerNode = FieldVariableType::element_type::BasisOnMesh::nDofsPerNode();
  
  // get all dofs that are associated with the current node
  std::array<dof_no_t,nDofsPerNode> dofGlobalNos;
  fieldVariable->mesh()->getNodeDofs(nodeGlobalNo, dofGlobalNos);
  
  // loop over components of the field variable
  for (int componentNo = 0; componentNo < fieldVariable->getNComponents(); componentNo++)
  {
    // get the values of the dofs for the current component
    std::array<double,nDofsPerNode> values;
    fieldVariable->template getValues<nDofsPerNode>(componentNo, dofGlobalNos, values);
    
    // output values
    for (double value : values)
    {
      stream << "  " << std::scientific << std::setprecision(17) << value << std::endl;
    }
  }
  
  // advance to next tuple element (next field variable)
  outputNodeValues<OutputFieldVariablesType, i+1>(fieldVariables, stream, nodeGlobalNo);
}

};   // namespace

 
//! write exnode file to given stream
template<typename BasisOnMeshType, typename OutputFieldVariablesType>
void ExfileWriter<BasisOnMeshType,OutputFieldVariablesType>::
outputExelem(std::ostream &stream, OutputFieldVariablesType fieldVariables)
{
  const int D = BasisOnMeshType::dim();
  stream << " Group name: Region" << std::endl
    << " Shape. Dimension=" << D << ", " << StringUtility::multiply<D>("line") << std::endl
    << " #Scale factor sets=0" << std::endl;
  
  const int nNodesPerElement = BasisOnMeshType::nNodesPerElement();   
  const element_no_t nElements = std::get<0>(fieldVariables)->mesh()->nElements();
  
  stream << " #Nodes=" << nNodesPerElement << std::endl
    << " #Fields=" << std::tuple_size<OutputFieldVariablesType>::value << std::endl;
     
  // loop over field variables and output headers
  ExfileLoopOverTuple::outputHeaderExelem<OutputFieldVariablesType>(fieldVariables, stream, 0);
  //int fieldVariableNo = 0;     // a number that runs over the field variables
  //for (auto &fieldVariable : fieldVariables)
  //{
//    fieldVariable->outputHeaderExelem(stream, 0, fieldVariableNo++);
  //}
      
  // loop over elements and output element-node information
/* example output
 Element:            1 0 0
   Nodes:
           1           2           4           5
*/           
  for(element_no_t elementGlobalNo = 0; elementGlobalNo < nElements; elementGlobalNo++)
  {
    stream << " Element:            " << elementGlobalNo+1 << " 0 0" << std::endl 
      << "   Nodes:" << std::endl;
     
    std::array<dof_no_t,BasisOnMeshType::nNodesPerElement()> elementNodes = std::get<0>(fieldVariables)->mesh()->getElementNodeNos(elementGlobalNo);
    StringUtility::outputValuesBlock(stream, elementNodes.begin(), elementNodes.end());
  }
  
}

//! write exnode file to given stream
template<typename BasisOnMeshType, typename OutputFieldVariablesType>
void ExfileWriter<BasisOnMeshType,OutputFieldVariablesType>::
outputExnode(std::ostream &stream, OutputFieldVariablesType fieldVariables)
{
  stream << " Group name: Region" << std::endl
    << " #Fields=" << std::tuple_size<OutputFieldVariablesType>::value << std::endl;
    
  // loop over field variables and output headers
  int valueIndex = 0;
  ExfileLoopOverTuple::outputHeaderExnode<OutputFieldVariablesType>(fieldVariables, stream, 0, valueIndex);
    
  //int fieldVariableNo = 0;     // a number that runs over the field variables
  //for (auto &fieldVariable : fieldVariables)
  //{
//    //(std::ostream &file, node_no_t currentNodeGlobalNo, int &valueIndex, int fieldVariableNo)
    //fieldVariable->outputHeaderExnode(stream, 0, valueIndex, fieldVariableNo++);
  //}

  // loop over nodes and output values
  const int nNodes = std::get<0>(fieldVariables)->mesh()->nNodes();
  for (node_no_t nodeGlobalNo = 0; nodeGlobalNo < nNodes; nodeGlobalNo++)
  {
    stream << " Node: " << nodeGlobalNo << std::endl;
    
    ExfileLoopOverTuple::outputNodeValues<OutputFieldVariablesType>(fieldVariables, stream, nodeGlobalNo);
    /*
    // loop over field variables 
    for (auto &fieldVariableBase : fieldVariables)
    {
      // get all values of the element for the field variable
      const int nComponents = fieldVariableBase->getNComponents();
      std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,nComponents>> fieldVariable
       = std::static_pointer_cast<FieldVariable::FieldVariable<BasisOnMeshType,nComponents>>(fieldVariableBase);
      
      std::array<std::array<double,nComponents>,BasisOnMeshType::nDofsPerElement()> elementValues;
      
      fieldVariable->template getElementValues<nComponents>(nodeGlobalNo, elementValues);
      
      // output values 
      for (auto &componentValue : elementValues)
      {
        for (double value : componentValue)
        {
          stream << std::scientific << std::setprecision(17) << value << std::endl;
        }
      }
    }*/
  }
/*
 Node:            1
  0.0000000000000000E+00
  0.0000000000000000E+00
  0.0000000000000000E+00
  -1.1100000000000001E+01
*/
}

  
};  //namespace