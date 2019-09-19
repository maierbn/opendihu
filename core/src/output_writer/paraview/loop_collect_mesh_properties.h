#pragma once

#include "utility/type_utility.h"

#include <cstdlib>

#include "output_writer/paraview/poly_data_properties_for_mesh.h"

/** The functions in this file model a loop over the elements of a tuple, as it occurs as FieldVariablesForOutputWriterType in all data_management classes.
 *  (Because the types inside the tuple are static and fixed at compile-time, a simple for loop c not work here.)
 *  The two functions starting with loop recursively emulate the loop. One method is the break condition and does nothing, the other method does the work and calls the method without loop in the name.
 *  FieldVariablesForOutputWriterType is assumed to be of type std::tuple<...>> where the types can be (mixed) std::shared_ptr<FieldVariable> or std::vector<std::shared_ptr<FieldVariable>>.
 * 
 *  Collect all field variable names of field variable with > 1 components into vectors and all scalar field variable names into scalars.
 */

namespace OutputWriter
{

namespace ParaviewLoopOverTuple
{

using ::OutputWriter::PolyDataPropertiesForMesh;
 
 /** Static recursive loop from 0 to number of entries in the tuple
 *  Stopping criterion
 */
template<typename FieldVariablesForOutputWriterType, int i=0>
inline typename std::enable_if<i == std::tuple_size<FieldVariablesForOutputWriterType>::value, void>::type
loopCollectMeshProperties(const FieldVariablesForOutputWriterType &fieldVariables,
                          std::map<std::string,PolyDataPropertiesForMesh> &meshProperties
)
{}

 /** Static recursive loop from 0 to number of entries in the tuple
 * Loop body
 */
template<typename FieldVariablesForOutputWriterType, int i=0>
inline typename std::enable_if<i < std::tuple_size<FieldVariablesForOutputWriterType>::value, void>::type
loopCollectMeshProperties(const FieldVariablesForOutputWriterType &fieldVariables,
                          std::map<std::string,PolyDataPropertiesForMesh> &meshProperties);

/** Loop body for a vector element
 */
template<typename VectorType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<TypeUtility::isVector<VectorType>::value, bool>::type
collectMeshProperties(VectorType currentFieldVariableVector, const FieldVariablesForOutputWriterType &fieldVariables,
                           std::map<std::string,PolyDataPropertiesForMesh> &meshProperties);

/** Loop body for a tuple element
 */
template<typename VectorType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<TypeUtility::isTuple<VectorType>::value, bool>::type
collectMeshProperties(VectorType currentFieldVariableVector, const FieldVariablesForOutputWriterType &fieldVariables,
                           std::map<std::string,PolyDataPropertiesForMesh> &meshProperties);

 /**  Loop body for a pointer element
 */
template<typename CurrentFieldVariableType, typename FieldVariablesForOutputWriterType>
typename std::enable_if<!TypeUtility::isTuple<CurrentFieldVariableType>::value && !TypeUtility::isVector<CurrentFieldVariableType>::value, bool>::type
collectMeshProperties(CurrentFieldVariableType currentFieldVariable, const FieldVariablesForOutputWriterType &fieldVariables,
                           std::map<std::string,PolyDataPropertiesForMesh> &meshProperties);

}  // namespace ParaviewLoopOverTuple

}  // namespace OutputWriter

#include "output_writer/paraview/loop_collect_mesh_properties.tpp"
