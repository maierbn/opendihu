#pragma once

#include "utility/type_utility.h"

#include <cstdlib>

/** The functions in this file model a loop over the nodeents of a tuple, as it occurs as OutputPointDataFieldVariablesType in all data_management classes.
 *  (Because the types inside the tuple are static and fixed at compile-time, a simple for loop c not work here.)
 *  The two functions starting with loop recursively emulate the loop. One method is the break condition and does nothing, the other method does the work and calls the method without loop in the name.
 *  OutputPointDataFieldVariablesType is assumed to be of type std::tuple<...>> where the types can be (mixed) std::shared_ptr<FieldVariable> or std::vector<std::shared_ptr<FieldVariable>>.
 * 
 *  Call ParaviewWriter::outputPointData on the mesh with meshName. This outputPointDatas all field variables of the mesh to a paraview readable file.
 */

namespace OutputWriter
{

namespace ParaviewLoopOverTuple
{
 
 /** Static recursive loop from 0 to number of entries in the tuple
 *  Stopping criterion
 */
template<typename OutputPointDataFieldVariablesType, int i=0>
inline typename std::enable_if<i == std::tuple_size<OutputPointDataFieldVariablesType>::value, void>::type
loopOutputPointData(const OutputPointDataFieldVariablesType &fieldVariables, std::string meshName,
                    std::ofstream &file, bool binaryOutput, bool fixedFormat
)
{}

 /** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename OutputPointDataFieldVariablesType, int i=0>
inline typename std::enable_if<i < std::tuple_size<OutputPointDataFieldVariablesType>::value, void>::type
loopOutputPointData(const OutputPointDataFieldVariablesType &fieldVariables, std::string meshName, 
                    std::ofstream &file, bool binaryOutput, bool fixedFormat);

/** Loop body for a vector nodeent
 */
template<typename VectorType, typename OutputPointDataFieldVariablesType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
outputPointData(VectorType currentFieldVariableVector, const OutputPointDataFieldVariablesType &fieldVariables, std::string meshName, 
                std::ofstream &file, bool binaryOutput, bool fixedFormat);

 /**  Loop body for a pointer nodeent
 */
template<typename CurrentFieldVariableType, typename OutputPointDataFieldVariablesType>
typename std::enable_if<!TypeUtility::isVector<CurrentFieldVariableType>::value, bool>::type
outputPointData(CurrentFieldVariableType currentFieldVariable, const OutputPointDataFieldVariablesType &fieldVariables, std::string meshName, 
                std::ofstream &file, bool binaryOutput, bool fixedFormat);

};  //namespace ParaviewLoopOverTuple

};  //namespace OutputWriter

#include "output_writer/paraview/loop_output_point_data.tpp"