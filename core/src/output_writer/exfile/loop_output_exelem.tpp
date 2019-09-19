#include "output_writer/exfile/loop_output_exelem.h"

#include "output_writer/loop_count_n_field_variables_of_mesh.h"

#include <cstdlib>

namespace OutputWriter
{

namespace ExfileLoopOverTuple
{
 
 /** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename FieldVariablesForOutputWriterType, typename AllFieldVariablesForOutputWriterType, int i>
inline typename std::enable_if<i < std::tuple_size<FieldVariablesForOutputWriterType>::value, void>::type
loopOutputExelem(const FieldVariablesForOutputWriterType &fieldVariables, const AllFieldVariablesForOutputWriterType &allFieldVariables,
                 std::string meshName, std::ofstream &file, std::shared_ptr<Mesh::Mesh> &mesh
)
{
  // call what to do in the loop body
  if (outputExelem<typename std::tuple_element<i,FieldVariablesForOutputWriterType>::type, AllFieldVariablesForOutputWriterType>(
        std::get<i>(fieldVariables), allFieldVariables, meshName, file, mesh))
    return;
  
  // advance iteration to next tuple element
  loopOutputExelem<FieldVariablesForOutputWriterType, AllFieldVariablesForOutputWriterType, i+1>(fieldVariables, allFieldVariables, meshName, file, mesh);
}
 
// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value, bool>::type
outputExelem(CurrentFieldVariableType currentFieldVariable, const FieldVariablesForOutputWriterType &fieldVariables, std::string meshName, 
             std::ofstream &file, std::shared_ptr<Mesh::Mesh> &mesh)
{
  // if mesh name is the specified meshName
  if (currentFieldVariable->functionSpace()->meshName() == meshName)
  {
    // here we have the type of the mesh with meshName (which is typedef to FunctionSpace)
    typedef typename CurrentFieldVariableType::element_type::FunctionSpace FunctionSpace;
   
    int nFieldVariablesInMesh = 0;
    LoopOverTuple::loopCountNFieldVariablesOfMesh(fieldVariables, meshName, nFieldVariablesInMesh);
    // get mesh
    mesh = currentFieldVariable->functionSpace();
    
    // call exfile writer to output all field variables with the meshName
    ExfileWriter<FunctionSpace, FieldVariablesForOutputWriterType>::outputExelem(file, fieldVariables, meshName, currentFieldVariable->functionSpace(), nFieldVariablesInMesh);
   
    return true;  // break iteration
  }
  
  return false;  // do not break iteration 
}

// Elementent i is of vector type
template<typename VectorType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
outputExelem(VectorType currentFieldVariableVector, const FieldVariablesForOutputWriterType &fieldVariables, std::string meshName, 
             std::ofstream &file, std::shared_ptr<Mesh::Mesh> &mesh)
{
  for (auto& currentFieldVariable : currentFieldVariableVector)
  {
    // call function on all vector entries
    if (outputExelem<typename VectorType::value_type,FieldVariablesForOutputWriterType>(currentFieldVariable, fieldVariables, meshName, file, mesh))
      return true; // break iteration
  }
  return false;  // do not break iteration 
}

// element i is of tuple type
template<typename TupleType, typename AllFieldVariablesForOutputWriterType>
typename std::enable_if<TypeUtility::isTuple<TupleType>::value, bool>::type
outputExelem(TupleType currentFieldVariableTuple, const AllFieldVariablesForOutputWriterType &fieldVariables, std::string meshName, 
             std::ofstream &file, std::shared_ptr<Mesh::Mesh> &mesh)
{
  // call for tuple element
  loopOutputExelem<TupleType, AllFieldVariablesForOutputWriterType>(currentFieldVariableTuple, fieldVariables, meshName, file, mesh);
  
  return false;  // do not break iteration 
}

}  // namespace ExfileLoopOverTuple
}  // namespace OutputWriter
