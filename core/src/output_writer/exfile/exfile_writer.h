#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>
#include <memory>

#include "basis_on_mesh/05_basis_on_mesh.h"

namespace OutputWriter
{

/** Base class of ExfileWriter that writes exelem and exnode files of given field variables.
 *  OutputFieldVariablesType is a std::tuple<std::shared_ptr<>, std::shared_ptr<>, ...> of field variables that will be output.
 */
template<typename BasisOnMeshType, typename OutputFieldVariablesType>
class ExfileWriter
{
public:
 
  //! write exnode file to given stream
  static void outputExelem(std::ostream &stream, OutputFieldVariablesType fieldVariables);
  
  //! write exnode file to given stream
  static void outputExnode(std::ostream &stream, OutputFieldVariablesType fieldVariables);
};

/*
// specialization for RegularFixed
template<int D, typename BasisFunctionType, typename OutputFieldVariablesType>
class ExfileWriter<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>, OutputFieldVariablesType>
{
public:
  typedef BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType> BasisOnMeshType;
 
  //! write exnode file to given stream
  static void outputExelem(std::ostream &stream, OutputFieldVariablesType fieldVariables);
  
  //! write exnode file to given stream
  static void outputExnode(std::ostream &stream, OutputFieldVariablesType fieldVariables);
};


// specialization for StructuredDeformable 
template<int D, typename BasisFunctionType, typename OutputFieldVariablesType>
class ExfileWriter<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>, OutputFieldVariablesType>
{
public:
  typedef BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType> BasisOnMeshType;
  
  //! write exnode file to given stream
  static void outputExelem(std::ostream &stream, OutputFieldVariablesType fieldVariables);
  
  //! write exnode file to given stream
  static void outputExnode(std::ostream &stream, OutputFieldVariablesType fieldVariables);
};

*/
// specialization for UnstructuredDeformable
template<int D, typename BasisFunctionType, typename OutputFieldVariablesType>
class ExfileWriter<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>, OutputFieldVariablesType>
{
public:
  typedef BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType> BasisOnMeshType;
  
  //! write exnode file to given stream
  static void outputExelem(std::ostream &stream, OutputFieldVariablesType fieldVariables);
  
  //! write exnode file to given stream
  static void outputExnode(std::ostream &stream, OutputFieldVariablesType fieldVariables);
};

namespace ExfileLoopOverTuple
{

 /** Static recursive loop from 0 to number of entries in the tuple
 *  Stopping criterion
 */
template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i == std::tuple_size<OutputFieldVariablesType>::value, void>::type
outputHeaderExelem(const OutputFieldVariablesType &fieldVariables, std::ostream &stream, element_no_t currentElementGlobalNo)
{}

/** Static recursive loop from 0 to number of entries in the tuple
 *  Loop body
 */
template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i < std::tuple_size<OutputFieldVariablesType>::value, void>::type
outputHeaderExelem(const OutputFieldVariablesType &fieldVariables, std::ostream &stream, element_no_t currentElementGlobalNo);


 /** Static recursive loop from 0 to number of entries in the tuple
 *  Stopping criterion
 */
template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i == std::tuple_size<OutputFieldVariablesType>::value, void>::type
outputHeaderExnode(const OutputFieldVariablesType &fieldVariables, std::ostream &stream, node_no_t currentNodeGlobalNo, int &valueIndex)
{}

/** Static recursive loop from 0 to number of entries in the tuple
 *  Loop body
 */
template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i < std::tuple_size<OutputFieldVariablesType>::value, void>::type
outputHeaderExnode(const OutputFieldVariablesType &fieldVariables, std::ostream &stream, node_no_t currentNodeGlobalNo, int &valueIndex);

 /** Static recursive loop from 0 to number of entries in the tuple
 *  Stopping criterion
 */
template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i == std::tuple_size<OutputFieldVariablesType>::value, void>::type
outputNodeValues(const OutputFieldVariablesType &fieldVariables, std::ostream &stream, node_no_t nodeGlobalNo)
{}

/** Static recursive loop from 0 to number of entries in the tuple
 *  Loop body
 */
template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i < std::tuple_size<OutputFieldVariablesType>::value, void>::type
outputNodeValues(const OutputFieldVariablesType &fieldVariables, std::ostream &stream, node_no_t nodeGlobalNo);

 /** Static recursive loop from 0 to number of entries in the tuple
 *  Stopping criterion
 */
template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i == std::tuple_size<OutputFieldVariablesType>::value, void>::type
checkIfNewExelemHeaderNecessary(const OutputFieldVariablesType &fieldVariables, element_no_t currentElementGlobalNo, bool &newHeaderNecessary)
{}

/** Static recursive loop from 0 to number of entries in the tuple
 *  Loop body
 */
template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i < std::tuple_size<OutputFieldVariablesType>::value, void>::type
checkIfNewExelemHeaderNecessary(const OutputFieldVariablesType &fieldVariables, element_no_t currentElementGlobalNo, bool &newHeaderNecessary);

 /** Static recursive loop from 0 to number of entries in the tuple
 *  Stopping criterion
 */
template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i == std::tuple_size<OutputFieldVariablesType>::value, void>::type
checkIfNewExnodeHeaderNecessary(const OutputFieldVariablesType &fieldVariables, element_no_t currentNodeGlobalNo, bool &newHeaderNecessary)
{}

/** Static recursive loop from 0 to number of entries in the tuple
 *  Loop body
 */
template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i < std::tuple_size<OutputFieldVariablesType>::value, void>::type
checkIfNewExnodeHeaderNecessary(const OutputFieldVariablesType &fieldVariables, element_no_t currentNodeGlobalNo, bool &newHeaderNecessary);

 /** Static recursive loop from 0 to number of entries in the tuple
 *  Stopping criterion
 */
template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i == std::tuple_size<OutputFieldVariablesType>::value, void>::type
getValuesAtNode(const OutputFieldVariablesType &fieldVariables, element_no_t currentNodeGlobalNo, std::vector<double> &valuesAtNode)
{}

/** Static recursive loop from 0 to number of entries in the tuple
 *  Loop body
 */
template<typename OutputFieldVariablesType, int i=0>
inline typename std::enable_if<i < std::tuple_size<OutputFieldVariablesType>::value, void>::type
getValuesAtNode(const OutputFieldVariablesType &fieldVariables, element_no_t currentNodeGlobalNo, std::vector<double> &valuesAtNode);

};   // namespace

};  // namespace

#include "output_writer/exfile/exfile_writer_structured.tpp"
#include "output_writer/exfile/exfile_writer_unstructured_deformable.tpp"
