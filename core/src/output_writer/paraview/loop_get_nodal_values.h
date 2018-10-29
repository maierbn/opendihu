#pragma once

#include "utility/type_utility.h"

#include <cstdlib>

/** The functions in this file model a loop over the elements of a tuple, as it occurs as OutputFieldVariablesType in all data_management classes.
 *  (Because the types inside the tuple are static and fixed at compile-time, a simple for loop c not work here.)
 *  The two functions starting with loop recursively emulate the loop. One method is the break condition and does nothing, the other method does the work and calls the method without loop in the name.
 *  OutputFieldVariablesType is assumed to be of type std::tuple<...>> where the types can be (mixed) std::shared_ptr<FieldVariable> or std::vector<std::shared_ptr<FieldVariable>>.
 * 
 *  Get all nodal values of all field variables in the meshes given in meshNames. The output variable is a vector of the data for the field variables (values[fieldVariableNo][entryNo]).
 *  The components of the field variable are interleaved (x y z x y z ...) as needed in the paraview output files.
 */

namespace OutputWriter
{

namespace ParaviewLoopOverTuple
{
 
 /** Static recursive loop from 0 to number of entries in the tuple
 *  Stopping criterion
 */
template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i == std::tuple_size<OutputFieldVariablesType>::value, void>::type
loopGetNodalValues(const OutputFieldVariablesType &fieldVariables, std::set<std::string> meshNames,
                   std::vector<std::vector<double>> &values
)
{}

 /** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i < std::tuple_size<OutputFieldVariablesType>::value, void>::type
loopGetNodalValues(const OutputFieldVariablesType &fieldVariables, std::set<std::string> meshNames,
                   std::vector<std::vector<double>> &values);

/** Loop body for a vector element
 */
template<typename VectorType, typename OutputFieldVariablesType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
getNodalValues(VectorType currentFieldVariableVector, const OutputFieldVariablesType &fieldVariables, std::set<std::string> meshNames,
               std::vector<std::vector<double>> &values);

/** Loop body for a tuple element
 */
template<typename VectorType, typename OutputFieldVariablesType>
typename std::enable_if<TypeUtility::isTuple<VectorType>::value, bool>::type
getNodalValues(VectorType currentFieldVariableVector, const OutputFieldVariablesType &fieldVariables, std::set<std::string> meshNames,
               std::vector<std::vector<double>> &values);

 /**  Loop body for a pointer element
 */
template<typename CurrentFieldVariableType, typename OutputFieldVariablesType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value, bool>::type
getNodalValues(CurrentFieldVariableType currentFieldVariable, const OutputFieldVariablesType &fieldVariables, std::set<std::string> meshNames,
               std::vector<std::vector<double>> &values);

};  //namespace ParaviewLoopOverTuple

};  //namespace OutputWriter

#include "output_writer/paraview/loop_get_nodal_values.tpp"
