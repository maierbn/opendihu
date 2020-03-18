#include "output_writer/exfile/loop_output_header_exnode.h"

#include <cstdlib>
#include "field_variable/field_variable.h"

namespace OutputWriter
{

namespace ExfileLoopOverTuple
{
 
 /** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename FieldVariablesForOutputWriterType, int i>
inline typename std::enable_if<i < std::tuple_size<FieldVariablesForOutputWriterType>::value, void>::type
loopOutputHeaderExnode(const FieldVariablesForOutputWriterType &fieldVariables, int &fieldVariableIndex, std::string meshName, 
                                      std::ostream &stream, node_no_t currentNodeGlobalNo, int &valueIndex)
{
  // call what to do in the loop body
  if (outputHeaderExnode<typename std::tuple_element<i,FieldVariablesForOutputWriterType>::type>(
        std::get<i>(fieldVariables), fieldVariableIndex, meshName, stream, currentNodeGlobalNo, valueIndex))
    return;
  
  // advance iteration to next tuple element
  loopOutputHeaderExnode<FieldVariablesForOutputWriterType, i+1>(fieldVariables, fieldVariableIndex, meshName, stream, currentNodeGlobalNo, valueIndex);
}
 
// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value && !Mesh::isComposite<CurrentFieldVariableType>::value, bool>::type
outputHeaderExnode(CurrentFieldVariableType currentFieldVariable, int &fieldVariableIndex, std::string meshName, 
                   std::ostream &stream, node_no_t currentNodeGlobalNo, int &valueIndex)
{
  // if mesh name is not the specified meshName step over this field variable but do not exit the loop over field variables
  if (currentFieldVariable->functionSpace()->meshName() != meshName)
  {
    return false;  // do not break iteration
  }
  
  currentFieldVariable->outputHeaderExnode(stream, currentNodeGlobalNo, valueIndex, fieldVariableIndex);
  fieldVariableIndex++;
  
  return false;  // do not break iteration 
}

// element i is of vector type
template<typename VectorType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
outputHeaderExnode(VectorType currentFieldVariableVector, int &fieldVariableIndex, std::string meshName, 
                   std::ostream &stream, node_no_t currentNodeGlobalNo, int &valueIndex)
{
  for (auto& currentFieldVariable : currentFieldVariableVector)
  {
    // call function on all vector entries
    if (outputHeaderExnode<typename VectorType::value_type>(currentFieldVariable, fieldVariableIndex, meshName, stream, currentNodeGlobalNo, valueIndex))
      return true;
  }
  return false;  // do not break iteration 
}

// element i is of tuple type
template<typename TupleType>
typename std::enable_if<TypeUtility::isTuple<TupleType>::value, bool>::type
outputHeaderExnode(TupleType currentFieldVariableTuple, int &fieldVariableIndex, std::string meshName, 
                   std::ostream &stream, node_no_t currentNodeGlobalNo, int &valueIndex)
{
  // call for tuple element
  loopOutputHeaderExnode<TupleType>(currentFieldVariableTuple, fieldVariableIndex, meshName, stream, currentNodeGlobalNo, valueIndex);
  
  return false;  // do not break iteration 
}

// element i is a field variables with Mesh::CompositeOfDimension<D>
template<typename CurrentFieldVariableType>
typename std::enable_if<Mesh::isComposite<CurrentFieldVariableType>::value, bool>::type
outputHeaderExnode(CurrentFieldVariableType currentFieldVariable, int &fieldVariableIndex, std::string meshName,
                   std::ostream &stream, node_no_t currentNodeGlobalNo, int &valueIndex)
{
  const int D = CurrentFieldVariableType::element_type::FunctionSpace::dim();
  typedef typename CurrentFieldVariableType::element_type::FunctionSpace::BasisFunction BasisFunctionType;
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<D>, BasisFunctionType> SubFunctionSpaceType;
  const int nComponents = CurrentFieldVariableType::element_type::nComponents();

  typedef FieldVariable::FieldVariable<SubFunctionSpaceType, nComponents> SubFieldVariableType;

  std::vector<std::shared_ptr<SubFieldVariableType>> subFieldVariables;
  currentFieldVariable->getSubFieldVariables(subFieldVariables);

  for (auto& currentSubFieldVariable : subFieldVariables)
  {
    // call function on all vector entries
    if (outputHeaderExnode<std::shared_ptr<SubFieldVariableType>>(currentSubFieldVariable, fieldVariableIndex, meshName, stream, currentNodeGlobalNo, valueIndex))
      return true;
  }

  return false;  // do not break iteration
}
}  // namespace ExfileLoopOverTuple
}  // namespace OutputWriter
