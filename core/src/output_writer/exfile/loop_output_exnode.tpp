#include "output_writer/exfile/loop_output_exnode.h"

#include "output_writer/loop_count_n_field_variables_of_mesh.h"

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
loopOutputExnode(const OutputFieldVariablesType &fieldVariables, std::string meshName, 
                 std::ofstream &file
)
{
  // call what to do in the loop body
  if (outputExnode<typename std::tuple_element<i,OutputFieldVariablesType>::type, OutputFieldVariablesType>(
        std::get<i>(fieldVariables), fieldVariables, meshName, file))
    return;
  
  // advance iteration to next tuple element
  loopOutputExnode<OutputFieldVariablesType, i+1>(fieldVariables, meshName, file);
}
 
// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType, typename OutputFieldVariablesType>
typename std::enable_if<!TypeUtility::isVector<CurrentFieldVariableType>::value, bool>::type
outputExnode(CurrentFieldVariableType currentFieldVariable, const OutputFieldVariablesType &fieldVariables, std::string meshName, 
             std::ofstream &file)
{
  // if mesh name is the specified meshName
  if (currentFieldVariable->mesh()->meshName() == meshName)
  {
    // here we have the type of the mesh with meshName (which is typedef to BasisOnMesh)
    typedef typename CurrentFieldVariableType::element_type::BasisOnMesh BasisOnMesh;
   
    // count number of field variables for the particular mesh
    int nFieldVariablesInMesh = 0;
    LoopOverTuple::loopCountNFieldVariablesOfMesh(fieldVariables, meshName, nFieldVariablesInMesh);
    
    // call exfile writer to output all field variables with the meshName
    ExfileWriter<BasisOnMesh, OutputFieldVariablesType>::outputExnode(file, fieldVariables, meshName, currentFieldVariable->mesh(), nFieldVariablesInMesh);
   
    return true;  // break iteration
  }
  
  return false;  // do not break iteration 
}

// Elementent i is of vector type
template<typename VectorType, typename OutputFieldVariablesType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
outputExnode(VectorType currentFieldVariableVector, const OutputFieldVariablesType &fieldVariables, std::string meshName, 
             std::ofstream &file)
{
  for (auto& currentFieldVariable : currentFieldVariableVector)
  {
    // call function on all vector entries
    if (outputExnode<typename VectorType::value_type,OutputFieldVariablesType>(currentFieldVariable, fieldVariables, meshName, file))
      return true; // break iteration
  }
  return false;  // do not break iteration 
}

};  //namespace ExfileLoopOverTuple
};  //namespace OutputWriter

