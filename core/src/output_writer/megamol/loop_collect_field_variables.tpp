#include "output_writer/megamol/loop_collect_field_variables.h"

#include <cstdlib>

namespace OutputWriter
{

namespace MegaMolLoopOverTuple
{
 
/** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename FieldVariablesForOutputWriterType, typename FunctionSpaceType, int i>
inline typename std::enable_if<i < std::tuple_size<FieldVariablesForOutputWriterType>::value, void>::type
loopCollectFieldVariables(const FieldVariablesForOutputWriterType &fieldVariables, std::string meshName,
                          std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> &geometryField,
                          std::vector<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>> &scalarFieldVariables)
{
  // call what to do in the loop body
  if (collectFieldVariables<typename std::tuple_element<i,FieldVariablesForOutputWriterType>::type,FunctionSpaceType>(std::get<i>(fieldVariables), meshName, geometryField, scalarFieldVariables))
    return;
  
  // advance iteration to next tuple element
  loopCollectFieldVariables<FieldVariablesForOutputWriterType, FunctionSpaceType, i+1>(fieldVariables, meshName, geometryField, scalarFieldVariables);
}
 
// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType, typename FunctionSpaceType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value
                        && CurrentFieldVariableType::element_type::nComponents() == 1, bool>::type
collectFieldVariables(CurrentFieldVariableType currentFieldVariable, std::string meshName,
                      std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> &geometryField,
                      std::vector<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>> &scalarFieldVariables)
{
  // if mesh name is not the specified meshName step over this field variable but do not exit the loop over field variables
  if (currentFieldVariable->functionSpace()->meshName() != meshName)
  {
    return false;  // do not break iteration
  }

  geometryField = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,3>>(currentFieldVariable->functionSpace()->geometryField());

  // number of components of currentFieldVariable is 1
  scalarFieldVariables.push_back(currentFieldVariable);

  return false;  // do not break iteration
}

// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType, typename FunctionSpaceType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value
                        && CurrentFieldVariableType::element_type::nComponents() != 1, bool>::type
collectFieldVariables(CurrentFieldVariableType currentFieldVariable, std::string meshName,
                      std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> &geometryField,
                      std::vector<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>> &scalarFieldVariables)
{
  // if mesh name is not the specified meshName step over this field variable but do not exit the loop over field variables
  if (currentFieldVariable->functionSpace()->meshName() != meshName)
  {
    return false;  // do not break iteration
  }

  geometryField = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,3>>(currentFieldVariable->functionSpace()->geometryField());

  return false;  // do not break iteration
}


// element i is of tuple type
template<typename TupleType, typename FunctionSpaceType>
typename std::enable_if<TypeUtility::isTuple<TupleType>::value, bool>::type
collectFieldVariables(TupleType currentFieldVariableTuple, std::string meshName,
                      std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> &geometryField,
                      std::vector<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>> &scalarFieldVariables)
{
  // call for tuple element
  loopCollectFieldVariables<TupleType,FunctionSpaceType>(currentFieldVariableTuple, meshName, geometryField, scalarFieldVariables);
  
  return false;  // do not break iteration
}

// element i is of vector type
template<typename VectorType, typename FunctionSpaceType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
collectFieldVariables(VectorType currentFieldVariableVector, std::string meshName,
                      std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> &geometryField,
                      std::vector<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>> &scalarFieldVariables)
{
  for (auto& currentFieldVariable : currentFieldVariableVector)
  {
    // call function on all vector entries
    if (collectFieldVariables<typename VectorType::value_type>(
     currentFieldVariable, meshName, geometryField, scalarFieldVariables))
      return true; // break iteration
  }
  return false;  // do not break iteration
}

}  // namespace MegaMolLoopOverTuple
}  // namespace OutputWriter
