#include "output_writer/megamol/loop_output.h"

#include "output_writer/megamol/megamol_writer.h"

#include <cstdlib>

namespace OutputWriter
{

namespace MegaMolLoopOverTuple
{
 
 /** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename OutputFieldVariablesType, typename AllOutputFieldVariablesType, int i>
inline typename std::enable_if<i < std::tuple_size<OutputFieldVariablesType>::value, void>::type
loopOutput(const OutputFieldVariablesType &fieldVariables, const AllOutputFieldVariablesType &allFieldVariables,
           std::string meshName, PythonConfig specificSettings, std::shared_ptr<adios2::Engine> adiosWriter, std::shared_ptr<adios2::IO> adiosIo, BoundingBox &boundingBox
)
{
  // call what to do in the loop body
  if (output<typename std::tuple_element<i,OutputFieldVariablesType>::type, AllOutputFieldVariablesType>(
        std::get<i>(fieldVariables), allFieldVariables, meshName, specificSettings, adiosWriter, adiosIo, boundingBox))
    return;
  
  // advance iteration to next tuple element
  loopOutput<OutputFieldVariablesType, AllOutputFieldVariablesType, i+1>(fieldVariables, allFieldVariables, meshName, specificSettings, adiosWriter, adiosIo, boundingBox);
}
 
// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType, typename OutputFieldVariablesType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value, bool>::type
output(CurrentFieldVariableType currentFieldVariable, const OutputFieldVariablesType &fieldVariables, std::string meshName, 
       PythonConfig specificSettings, std::shared_ptr<adios2::Engine> adiosWriter, std::shared_ptr<adios2::IO> adiosIo, BoundingBox &boundingBox)
{
  // if mesh name is the specified meshName
  if (currentFieldVariable->functionSpace()->meshName() == meshName)
  {
    // here we have the type of the mesh with meshName (which is typedef to FunctionSpace)
    typedef typename CurrentFieldVariableType::element_type::FunctionSpace FunctionSpace;
   
    // call exfile writer to output all field variables with the meshName
    MegaMolWriter<FunctionSpace, OutputFieldVariablesType>::outputData(fieldVariables, meshName, currentFieldVariable->functionSpace(), specificSettings, adiosWriter, adiosIo, boundingBox);
   
    return true;  // break iteration
  }
  
  return false;  // do not break iteration 
}

// element i is of vector type
template<typename VectorType, typename OutputFieldVariablesType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
output(VectorType currentFieldVariableVector, const OutputFieldVariablesType &fieldVariables, std::string meshName, 
       PythonConfig specificSettings, std::shared_ptr<adios2::Engine> adiosWriter, std::shared_ptr<adios2::IO> adiosIo, BoundingBox &boundingBox)
{
  for (auto& currentFieldVariable : currentFieldVariableVector)
  {
    // call function on all vector entries
    if (output<typename VectorType::value_type,OutputFieldVariablesType>(currentFieldVariable, fieldVariables, meshName, specificSettings, adiosWriter, adiosIo, boundingBox))
      return true; // break iteration
  }
  return false;  // do not break iteration 
}

// element i is of tuple type
template<typename TupleType, typename AllOutputFieldVariablesType>
typename std::enable_if<TypeUtility::isTuple<TupleType>::value, bool>::type
output(TupleType currentFieldVariableTuple, const AllOutputFieldVariablesType &fieldVariables, std::string meshName, 
      PythonConfig specificSettings, std::shared_ptr<adios2::Engine> adiosWriter, std::shared_ptr<adios2::IO> adiosIo, BoundingBox &boundingBox)
{
  // call for tuple element
  loopOutput<TupleType, AllOutputFieldVariablesType>(currentFieldVariableTuple, fieldVariables, meshName, specificSettings, adiosWriter, adiosIo, boundingBox);
  
  return false;  // do not break iteration 
}

}  // namespace MegaMolLoopOverTuple
}  // namespace OutputWriter
