#include "output_writer/exfile/loop_check_if_new_exnode_header_necessary.h"

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
typename std::enable_if<i < std::tuple_size<FieldVariablesForOutputWriterType>::value, void>::type
loopCheckIfNewExnodeHeaderNecessary(const FieldVariablesForOutputWriterType &fieldVariables, std::string meshName,
                                    element_no_t currentNodeGlobalNo, bool &newHeaderNecessary)
{
  // call what to do in the loop body
  if (checkIfNewExnodeHeaderNecessary<typename std::tuple_element<i,FieldVariablesForOutputWriterType>::type>(
        std::get<i>(fieldVariables), meshName, currentNodeGlobalNo, newHeaderNecessary))
    return;
  
  // advance iteration to next tuple element
  loopCheckIfNewExnodeHeaderNecessary<FieldVariablesForOutputWriterType, i+1>(fieldVariables, meshName, currentNodeGlobalNo, newHeaderNecessary);
}
 
// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value && !Mesh::isComposite<CurrentFieldVariableType>::value, bool>::type
checkIfNewExnodeHeaderNecessary(CurrentFieldVariableType currentFieldVariable, std::string meshName,
                                element_no_t currentNodeGlobalNo, bool &newHeaderNecessary
)
{
  // if mesh name is not the specified meshName step over this field variable but do not exit the loop over field variables
  if (currentFieldVariable->functionSpace()->meshName() != meshName)
  {
    return false;  // do not break iteration
  }
  
  int previousNumberVersions = currentFieldVariable->nodeToDofMapping()->nVersions(currentNodeGlobalNo-1);
  int currentNumberVersions = currentFieldVariable->nodeToDofMapping()->nVersions(currentNodeGlobalNo);

  if (previousNumberVersions != currentNumberVersions)
  {
    newHeaderNecessary = true;
  
    // break iteration
    return true;
  }
  
  return false;
}

// element i is of vector type
template<typename VectorType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
checkIfNewExnodeHeaderNecessary(VectorType currentFieldVariableVector, std::string meshName,
                                element_no_t currentNodeGlobalNo, bool &newHeaderNecessary)
{
  for (auto& currentFieldVariable : currentFieldVariableVector)
  {
    // call function on all vector entries
    if (checkIfNewExnodeHeaderNecessary<typename VectorType::value_type>(currentFieldVariable, meshName, currentNodeGlobalNo, newHeaderNecessary))
      return true;
  }
  return false;  // do not break iteration
}

// element i is of tuple type
template<typename TupleType>
typename std::enable_if<TypeUtility::isTuple<TupleType>::value, bool>::type
checkIfNewExnodeHeaderNecessary(TupleType currentFieldVariableTuple, std::string meshName,
                                element_no_t currentNodeGlobalNo, bool &newHeaderNecessary)
{
  // call for tuple element
  loopCheckIfNewExnodeHeaderNecessary<TupleType>(currentFieldVariableTuple, meshName, currentNodeGlobalNo, newHeaderNecessary);
  
  return false;  // do not break iteration
}

// element i is a field variables with Mesh::CompositeOfDimension<D>
template<typename CurrentFieldVariableType>
typename std::enable_if<Mesh::isComposite<CurrentFieldVariableType>::value, bool>::type
checkIfNewExnodeHeaderNecessary(CurrentFieldVariableType currentFieldVariable, std::string meshName,
                                element_no_t currentNodeGlobalNo, bool &newHeaderNecessary)
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
    if (checkIfNewExnodeHeaderNecessary<std::shared_ptr<SubFieldVariableType>>(currentSubFieldVariable, meshName, currentNodeGlobalNo, newHeaderNecessary))
      return true;
  }

  return false;  // do not break iteration
}
}  // namespace ExfileLoopOverTuple
}  // namespace OutputWriter
