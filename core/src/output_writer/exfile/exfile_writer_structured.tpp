#include "output_writer/exfile/exfile_writer.h"

#include <iomanip>

namespace OutputWriter
{

//! write exnode file to given stream
template<typename BasisOnMeshType, typename OutputFieldVariablesType>
void ExfileWriter<BasisOnMeshType,OutputFieldVariablesType>::
outputExelem(std::ostream &stream, OutputFieldVariablesType fieldVariables, std::string meshName, std::shared_ptr<BasisOnMeshType> mesh, int nFieldVariablesOfMesh)
{
  const int D = BasisOnMeshType::dim();
  stream << " Group name: " << meshName << std::endl
    << " Shape. Dimension=" << D << ", " << StringUtility::multiply<D>("line") << std::endl
    << " #Scale factor sets=0" << std::endl;

  const int nNodesPerElement = BasisOnMeshType::nNodesPerElement();
  const element_no_t nElements = mesh->nElementsLocal();

  stream << " #Nodes=" << nNodesPerElement << std::endl
    << " #Fields=" << nFieldVariablesOfMesh << std::endl;

  // loop over field variables and output headers
  int fieldVariableIndex = 0;  // counter over field variables
  ExfileLoopOverTuple::loopOutputHeaderExelem<OutputFieldVariablesType>(fieldVariables, fieldVariableIndex, meshName, stream, 0);
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

    std::array<dof_no_t,BasisOnMeshType::nNodesPerElement()> elementNodes = mesh->getElementNodeNos(elementGlobalNo);
    StringUtility::outputValuesBlockAdd1(stream, elementNodes.begin(), elementNodes.end());
  }

}

//! write exnode file to given stream
template<typename BasisOnMeshType, typename OutputFieldVariablesType>
void ExfileWriter<BasisOnMeshType,OutputFieldVariablesType>::
outputExnode(std::ostream &stream, OutputFieldVariablesType fieldVariables, std::string meshName, 
             std::shared_ptr<BasisOnMeshType> mesh, int nFieldVariablesOfMesh)
{
  VLOG(2) << "ExfileWriter<Structured>::outputExnode, meshName: " << meshName << ",nFieldVariablesOfMesh:" << nFieldVariablesOfMesh;
    
  stream << " Group name: " << meshName << std::endl
    << " #Fields=" << nFieldVariablesOfMesh << std::endl;

  // loop over field variables and output headers
  int valueIndex = 0;
  int fieldVariableIndex = 0;  // counter over field variables
  ExfileLoopOverTuple::loopOutputHeaderExnode<OutputFieldVariablesType>(fieldVariables, fieldVariableIndex, meshName, stream, 0, valueIndex);

  //int fieldVariableNo = 0;     // a number that runs over the field variables
  //for (auto &fieldVariable : fieldVariables)
  //{
//    //(std::ostream &file, node_no_t currentNodeGlobalNo, int &valueIndex, int fieldVariableNo)
    //fieldVariable->outputHeaderExnode(stream, 0, valueIndex, fieldVariableNo++);
  //}

  // loop over nodes and output values
  const int nNodes = mesh->nNodesGlobal();
  for (node_no_t nodeGlobalNo = 0; nodeGlobalNo < nNodes; nodeGlobalNo++)
  {
    stream << " Node: " << nodeGlobalNo+1 << std::endl;

    ExfileLoopOverTuple::loopOutputNodeValues<OutputFieldVariablesType>(fieldVariables, meshName, stream, nodeGlobalNo);
    /*
    // loop over field variables
    for (auto &fieldVariableBase : fieldVariables)
    {
      // get all values of the element for the field variable
      const int nComponents = fieldVariableBase->nComponents();
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