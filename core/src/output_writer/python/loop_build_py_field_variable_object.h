#pragma once

#include "utility/type_utility.h"

#include <cstdlib>

/** The functions in this file model a loop over the elements of a tuple, as it occurs as OutputFieldVariablesType in all data_management classes.
 *  (Because the types inside the tuple are static and fixed at compile-time, a simple for loop c not work here.)
 *  The two functions starting with loop recursively emulate the loop. One method is the break condition and does nothing, the other method does the work and calls the method without loop in the name.
 *  OutputFieldVariablesType is assumed to be of type std::tuple<...>> where the types can be (mixed) std::shared_ptr<FieldVariable> or std::vector<std::shared_ptr<FieldVariable>>.
 */

namespace OutputWriter
{

namespace PythonLoopOverTuple
{
 
 /** Static recursive loop from 0 to number of entries in the tuple
 *  Stopping criterion
 */
template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i == std::tuple_size<OutputFieldVariablesType>::value, void>::type
loopBuildPyFieldVariableObject(const OutputFieldVariablesType &fieldVariables, int &fieldVariableIndex, std::string meshName, 
                               PyObject *pyData, bool onlyNodalValues, std::shared_ptr<Mesh::Mesh> &mesh)
{}

 /** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i < std::tuple_size<OutputFieldVariablesType>::value, void>::type
loopBuildPyFieldVariableObject(const OutputFieldVariablesType &fieldVariables, int &fieldVariableIndex, std::string meshName, 
                               PyObject *pyData, bool onlyNodalValues, std::shared_ptr<Mesh::Mesh> &mesh);

/** Loop body for a vector element
 */
template<typename VectorType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
buildPyFieldVariableObject(VectorType currentFieldVariableVector, int &fieldVariableIndex, std::string meshName, 
                           PyObject *pyData, bool onlyNodalValues, std::shared_ptr<Mesh::Mesh> &mesh);

/** Loop body for a tuple element
 */
template<typename VectorType>
typename std::enable_if<TypeUtility::isTuple<VectorType>::value, bool>::type
buildPyFieldVariableObject(VectorType currentFieldVariableVector, int &fieldVariableIndex, std::string meshName, 
                           PyObject *pyData, bool onlyNodalValues, std::shared_ptr<Mesh::Mesh> &mesh);

/**  Loop body for a pointer element
 */
template<typename CurrentFieldVariableType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value, bool>::type
buildPyFieldVariableObject(CurrentFieldVariableType currentFieldVariable, int &fieldVariableIndex, std::string meshName, 
                           PyObject *pyData, bool onlyNodalValues, std::shared_ptr<Mesh::Mesh> &mesh);

}  // namespace ExfileLoopOverTuple

}  // namespace OutputWriter

#include "output_writer/python/loop_build_py_field_variable_object.tpp"
