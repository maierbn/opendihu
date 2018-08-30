#include "output_writer/loop_collect_mesh_names.h"

#include <cstdlib>

namespace OutputWriter
{

namespace LoopOverTuple
{
 
 /** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i < std::tuple_size<OutputFieldVariablesType>::value, void>::type
loopCollectMeshNames(const OutputFieldVariablesType &fieldVariables, std::set<std::string> &meshNames)
{
  // call what to do in the loop body
  if (collectMeshNames<typename std::tuple_element<i,OutputFieldVariablesType>::type>(
       std::get<i>(fieldVariables), meshNames))
    return;
  
  // advance iteration to next tuple element
  loopCollectMeshNames<OutputFieldVariablesType, i+1>(fieldVariables, meshNames);
}
 
// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value, bool>::type
collectMeshNames(CurrentFieldVariableType currentFieldVariable, std::set<std::string> &meshNames)
{
  // get mesh name and insert into meshNames
  std::string meshName = currentFieldVariable->functionSpace()->meshName();
  meshNames.insert(meshName);
  
  return false;  // do not break iteration
}

// element i is of vector type
template<typename VectorType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
collectMeshNames(VectorType currentFieldVariableVector, std::set<std::string> &meshNames)
{
  for (auto& currentFieldVariable : currentFieldVariableVector)
  {
    // call function on all vector entries
    if (collectMeshNames<typename VectorType::value_type>(currentFieldVariable, meshNames))
      return true;
  }
  
  return false;  // do not break iteration
}

// element i is of tuple type
template<typename TupleType>
typename std::enable_if<TypeUtility::isTuple<TupleType>::value, bool>::type
collectMeshNames(TupleType currentFieldVariableTuple, std::set<std::string> &meshNames)
{
  // call for tuple element
  loopCollectMeshNames<TupleType>(currentFieldVariableTuple, meshNames);
 
  return false;  // do not break iteration
}

};  //namespace ExfileLoopOverTuple
};  //namespace OutputWriter
