#include "output_writer/paraview/loop_collect_field_variables_names.h"

#include "output_writer/paraview/paraview_writer.h"

#include <cstdlib>

namespace OutputWriter
{

namespace ParaviewLoopOverTuple
{
 
 /** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename FieldVariablesForOutputWriterType, int i>
inline typename std::enable_if<i < std::tuple_size<FieldVariablesForOutputWriterType>::value, void>::type
loopCollectFieldVariablesNames(const FieldVariablesForOutputWriterType &fieldVariables, std::string meshName, 
                               std::vector<std::string> &scalars, std::vector<std::string> &vectors
)
{
  // call what to do in the loop body
  if (collectFieldVariablesNames<typename std::tuple_element<i,FieldVariablesForOutputWriterType>::type, FieldVariablesForOutputWriterType>(
        std::get<i>(fieldVariables), fieldVariables, meshName, scalars, vectors))
    return;
  
  // advance iteration to next tuple element
  loopCollectFieldVariablesNames<FieldVariablesForOutputWriterType, i+1>(fieldVariables, meshName, scalars, vectors);
}
 
// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value, bool>::type
collectFieldVariablesNames(CurrentFieldVariableType currentFieldVariable, const FieldVariablesForOutputWriterType &fieldVariables, std::string meshName, 
                           std::vector<std::string> &scalars, std::vector<std::string> &vectors)
{
  // if mesh name is the specified meshName
  if (currentFieldVariable->functionSpace()->meshName() == meshName && !currentFieldVariable->isGeometryField())
  {
    if (currentFieldVariable->nComponents() == 1)
    {
      scalars.push_back(currentFieldVariable->name());
    }
    else 
    {
      vectors.push_back(currentFieldVariable->name());
    }
    
    return true;  // break iteration
  }
  
  return false;  // do not break iteration 
}

// element i is of vector type
template<typename VectorType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
collectFieldVariablesNames(VectorType currentFieldVariableVector, const FieldVariablesForOutputWriterType &fieldVariables, std::string meshName, 
                           std::vector<std::string> &scalars, std::vector<std::string> &vectors)
{
  for (auto& currentFieldVariable : currentFieldVariableVector)
  {
    // call function on all vector entries
    if (collectFieldVariablesNames<typename VectorType::value_type,FieldVariablesForOutputWriterType>(currentFieldVariable, fieldVariables, meshName, scalars, vectors))
      return true; // break iteration
  }
  return false;  // do not break iteration 
}

// element i is of tuple type
template<typename TupleType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<TypeUtility::isTuple<TupleType>::value, bool>::type
collectFieldVariablesNames(TupleType currentFieldVariableTuple, const FieldVariablesForOutputWriterType &fieldVariables, std::string meshName, 
                           std::vector<std::string> &scalars, std::vector<std::string> &vectors)
{
  // call for tuple element
  loopCollectFieldVariablesNames<TupleType>(currentFieldVariableTuple, meshName, scalars, vectors);
 
  return false;  // do not break iteration 
}

}  // namespace ParaviewLoopOverTuple
}  // namespace OutputWriter
