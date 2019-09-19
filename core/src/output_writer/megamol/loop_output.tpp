#include "output_writer/megamol/loop_output.h"

#include "output_writer/megamol/megamol_writer.h"

#include <cstdlib>

#ifdef HAVE_ADIOS
namespace OutputWriter
{

namespace MegaMolLoopOverTuple
{
 
 /** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename FieldVariablesForOutputWriterType, typename AllFieldVariablesForOutputWriterType, int i>
inline typename std::enable_if<i < std::tuple_size<FieldVariablesForOutputWriterType>::value, void>::type
loopOutput(const FieldVariablesForOutputWriterType &fieldVariables, const AllFieldVariablesForOutputWriterType &allFieldVariables,
           std::string meshName, PythonConfig specificSettings, MegaMolWriterContext &megaMolWriterContext
)
{
  // call what to do in the loop body
  if (output<typename std::tuple_element<i,FieldVariablesForOutputWriterType>::type, AllFieldVariablesForOutputWriterType>(
        std::get<i>(fieldVariables), allFieldVariables, meshName, specificSettings, megaMolWriterContext))
    return;
  
  // advance iteration to next tuple element
  loopOutput<FieldVariablesForOutputWriterType, AllFieldVariablesForOutputWriterType, i+1>(fieldVariables, allFieldVariables, meshName,
                                                                         specificSettings, megaMolWriterContext);
}
 
// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value, bool>::type
output(CurrentFieldVariableType currentFieldVariable, const FieldVariablesForOutputWriterType &fieldVariables, std::string meshName, 
       PythonConfig specificSettings, MegaMolWriterContext &megaMolWriterContext)
{
  // if mesh name is the specified meshName
  if (currentFieldVariable->functionSpace()->meshName() == meshName)
  {
    // here we have the type of the mesh with meshName (which is typedef to FunctionSpace)
    typedef typename CurrentFieldVariableType::element_type::FunctionSpace FunctionSpace;
/*
    if (nNodesGlobalAllMeshes.find(meshName) == nNodesGlobalAllMeshes.end())
    {
      nNodesGlobalAllMeshes[meshName] = 0;
    }

    nNodesGlobalAllMeshes[meshName] = std::max(nNodesGlobalAllMeshes[meshName], currentFieldVariable->nDofsGlobal());
  */
    // call megamol writer to output all field variables with the meshName
    MegaMolWriter<FunctionSpace, FieldVariablesForOutputWriterType>::outputData(fieldVariables, meshName, currentFieldVariable->functionSpace(),
                                                                       specificSettings, megaMolWriterContext);
   
    return true;  // break iteration
  }
  
  return false;  // do not break iteration 
}

// element i is of vector type
template<typename VectorType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
output(VectorType currentFieldVariableVector, const FieldVariablesForOutputWriterType &fieldVariables, std::string meshName, 
       PythonConfig specificSettings, MegaMolWriterContext &megaMolWriterContext)
{
  for (auto& currentFieldVariable : currentFieldVariableVector)
  {
    // call function on all vector entries
    if (output<typename VectorType::value_type,FieldVariablesForOutputWriterType>(
      currentFieldVariable, fieldVariables, meshName, specificSettings, megaMolWriterContext)
    )
      return true; // break iteration
  }
  return false;  // do not break iteration 
}

// element i is of tuple type
template<typename TupleType, typename AllFieldVariablesForOutputWriterType>
typename std::enable_if<TypeUtility::isTuple<TupleType>::value, bool>::type
output(TupleType currentFieldVariableTuple, const AllFieldVariablesForOutputWriterType &fieldVariables, std::string meshName, 
      PythonConfig specificSettings, MegaMolWriterContext &megaMolWriterContext)
{
  // call for tuple element
  loopOutput<TupleType, AllFieldVariablesForOutputWriterType>(currentFieldVariableTuple, fieldVariables, meshName, specificSettings, megaMolWriterContext);
  
  return false;  // do not break iteration 
}

}  // namespace MegaMolLoopOverTuple
}  // namespace OutputWriter
#endif
