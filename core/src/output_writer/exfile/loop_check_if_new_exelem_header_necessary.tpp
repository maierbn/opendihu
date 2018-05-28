#include "output_writer/exfile/loop_check_if_new_exelem_header_necessary.h"

#include <cstdlib>

namespace OutputWriter
{

namespace ExfileLoopOverTuple
{
 
 /** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i < std::tuple_size<OutputFieldVariablesType>::value, void>::type
loopCheckIfNewExelemHeaderNecessary(const OutputFieldVariablesType &fieldVariables, std::string meshName, element_no_t currentFieldVariableGlobalNo, bool &newHeaderNecessary)
{
  // call what to do in the loop body
  if (checkIfNewExelemHeaderNecessary<typename std::tuple_element<i,OutputFieldVariablesType>::type>(std::get<i>(fieldVariables), meshName, currentFieldVariableGlobalNo, newHeaderNecessary))
    return;
  
  // advance iteration to next tuple element
  loopCheckIfNewExelemHeaderNecessary<OutputFieldVariablesType, i+1>(fieldVariables, meshName, currentFieldVariableGlobalNo, newHeaderNecessary);
}
 
// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType>
typename std::enable_if<!TypeUtility::isVector<CurrentFieldVariableType>::value, bool>::type
checkIfNewExelemHeaderNecessary(CurrentFieldVariableType currentFieldVariable, std::string meshName,
                                element_no_t currentFieldVariableGlobalNo, bool &newHeaderNecessary)
{
  // if mesh name is not the specified meshName step over this field variable but do not exit the loop over field variables
  if (currentFieldVariable->mesh()->meshName() != meshName)
  {
    return false;  // do not break iteration
  }
  
  VLOG(2) << "check if field variable " << currentFieldVariable->name() << " has the same exfileRepr for elements "
    << currentFieldVariableGlobalNo-1 << " and " << currentFieldVariableGlobalNo;

  if (!currentFieldVariable->haveSameExfileRepresentation(currentFieldVariableGlobalNo-1, currentFieldVariableGlobalNo))
  {
    newHeaderNecessary = true;
    
    // break iteration
    return true;
  }
  return false;   // do not break iteration
}

// element i is of vector type
template<typename VectorType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
checkIfNewExelemHeaderNecessary(VectorType currentFieldVariableVector, std::string meshName, 
                                element_no_t currentFieldVariableGlobalNo, bool &newHeaderNecessary)
{
  for (auto& currentFieldVariable : currentFieldVariableVector)
  {
    // call function on all vector entries
    if (checkIfNewExelemHeaderNecessary<typename VectorType::value_type>(currentFieldVariable, meshName, currentFieldVariableGlobalNo, newHeaderNecessary))
      return true;
  }
  
  return false;  // do not break iteration
}

};  //namespace ExfileLoopOverTuple
};  //namespace OutputWriter
