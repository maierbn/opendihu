#include "output_writer/paraview/loop_output_point_data.h"

#include <cstdlib>

#include "output_writer/paraview/paraview.h"

namespace OutputWriter
{

// forward declaration
struct Paraview;
 
namespace ParaviewLoopOverTuple
{
 
 /** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename OutputPointDataFieldVariablesType, int i=0>
inline typename std::enable_if<i < std::tuple_size<OutputPointDataFieldVariablesType>::value, void>::type
loopOutputPointData(const OutputPointDataFieldVariablesType &fieldVariables, std::string meshName, 
                    std::ofstream &file, bool binaryOutput, bool fixedFormat
)
{
  // call what to do in the loop body
  if (outputPointData<typename std::tuple_element<i,OutputPointDataFieldVariablesType>::type, OutputPointDataFieldVariablesType>(
        std::get<i>(fieldVariables), fieldVariables, meshName, file, binaryOutput, fixedFormat))
    return;
  
  // advance iteration to next tuple element
  loopOutputPointData<OutputPointDataFieldVariablesType, i+1>(fieldVariables, meshName, file, binaryOutput, fixedFormat);
}
 
// current element is of pointer type (not vector)
template<typename CurrentFieldVariableType, typename OutputPointDataFieldVariablesType>
typename std::enable_if<!TypeUtility::isVector<CurrentFieldVariableType>::value, bool>::type
outputPointData(CurrentFieldVariableType currentFieldVariable, const OutputPointDataFieldVariablesType &fieldVariables, std::string meshName, 
                std::ofstream &file, bool binaryOutput, bool fixedFormat)
{
  // if mesh name is the specified meshName
  if (currentFieldVariable->mesh()->meshName() == meshName && !currentFieldVariable->isGeometryField())
  {
    Paraview::writeParaviewFieldVariable<typename CurrentFieldVariableType::element_type>(*currentFieldVariable, file, binaryOutput, fixedFormat);
  }
  
  return false;  // do not break iteration 
}

// element i is of vector type
template<typename VectorType, typename OutputPointDataFieldVariablesType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
outputPointData(VectorType currentFieldVariableVector, const OutputPointDataFieldVariablesType &fieldVariables, std::string meshName, 
                std::ofstream &file, bool binaryOutput, bool fixedFormat)
{
  for (auto& currentFieldVariable : currentFieldVariableVector)
  {
    // call function on all vector entries
    if (outputPointData<typename VectorType::value_type,OutputPointDataFieldVariablesType>(currentFieldVariable, fieldVariables, meshName, file, binaryOutput, fixedFormat))
      return true; // break iteration
  }
  return false;  // do not break iteration 
}

};  //namespace ParaviewLoopOverTuple
};  //namespace OutputWriter
