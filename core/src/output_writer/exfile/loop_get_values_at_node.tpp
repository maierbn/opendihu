#include "output_writer/exfile/loop_get_values_at_node.h"

#include <cstdlib>

namespace OutputWriter
{

namespace ExfileLoopOverTuple
{
 
 /** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename FieldVariablesForOutputWriterType, int i>
inline typename std::enable_if<i < std::tuple_size<FieldVariablesForOutputWriterType>::value, void>::type
loopGetValuesAtNode(const FieldVariablesForOutputWriterType &fieldVariables, std::string meshName, 
                     element_no_t currentNodeGlobalNo, std::vector<double> &valuesAtNode)
{
  // call what to do in the loop body
  if (getValuesAtNode<typename std::tuple_element<i,FieldVariablesForOutputWriterType>::type>(std::get<i>(fieldVariables), meshName, currentNodeGlobalNo, valuesAtNode))
    return;
  
  // advance iteration to next tuple element
  loopGetValuesAtNode<FieldVariablesForOutputWriterType, i+1>(fieldVariables, meshName, currentNodeGlobalNo, valuesAtNode);
}
 
// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value, bool>::type
getValuesAtNode(CurrentFieldVariableType currentFieldVariable, std::string meshName, 
                element_no_t currentNodeGlobalNo, std::vector<double> &valuesAtNode)
{
  // if mesh name is not the specified meshName step over this field variable but do not exit the loop over field variables
  if (currentFieldVariable->functionSpace()->meshName() != meshName)
  {
    return false;  // do not break iteration
  }
 
  VLOG(1) << "  fieldVariable \"" << currentFieldVariable->name() << "\", "
    << " nodeToDofMapping: " << currentFieldVariable->nodeToDofMapping();
  VLOG(1) << "  meshName: " << meshName << ", get dofsAtNode for currentNodeGlobalNo: " << currentNodeGlobalNo << ", n components: " << currentFieldVariable->nComponents();

  std::vector<dof_no_t> &dofsAtNode = currentFieldVariable->nodeToDofMapping()->getNodeDofs(currentNodeGlobalNo);

  // loop over components
  for (int componentNo = 0; componentNo < currentFieldVariable->nComponents(); componentNo++)
  {
    currentFieldVariable->getValues(componentNo, dofsAtNode, valuesAtNode);
    VLOG(1) << "   Component " << componentNo;
    VLOG(1) << "    dofsAtNode: " << dofsAtNode; 
    VLOG(1) << "    valuesAtNode: " << valuesAtNode;
  }
  
  return false;  // do not break iteration
}


// element i is of tuple type
template<typename TupleType>
typename std::enable_if<TypeUtility::isTuple<TupleType>::value, bool>::type
getValuesAtNode(TupleType currentFieldVariableTuple, std::string meshName, 
                element_no_t currentNodeGlobalNo, std::vector<double> &valuesAtNode)
{
  // call for tuple element
  loopGetValuesAtNode<TupleType>(currentFieldVariableTuple, meshName, currentNodeGlobalNo, valuesAtNode);
  
  return false;  // do not break iteration
}

// element i is of vector type
template<typename VectorType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
getValuesAtNode(VectorType currentFieldVariableVector, std::string meshName, 
                element_no_t currentNodeGlobalNo, std::vector<double> &valuesAtNode)
{
  for (auto& currentFieldVariable : currentFieldVariableVector)
  {
    // call function on all vector entries
    if (getValuesAtNode<typename VectorType::value_type>(
     currentFieldVariable, meshName, currentNodeGlobalNo, valuesAtNode))
      return true; // break iteration
  }
  return false;  // do not break iteration
}

}  // namespace ExfileLoopOverTuple
}  // namespace OutputWriter
