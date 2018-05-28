#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>
#include <memory>

#include "basis_on_mesh/basis_on_mesh.h"
#include "output_writer/exfile/loop_check_if_new_exelem_header_necessary.h"
#include "output_writer/exfile/loop_check_if_new_exnode_header_necessary.h"
#include "output_writer/exfile/loop_get_values_at_node.h"
#include "output_writer/exfile/loop_output_header_exelem.h"
#include "output_writer/exfile/loop_output_header_exnode.h"
#include "output_writer/exfile/loop_output_node_values.h"

namespace OutputWriter
{

/** Base class of ExfileWriter that writes exelem and exnode files of given field variables.
 *  OutputFieldVariablesType is a std::tuple<std::shared_ptr<>, std::shared_ptr<>, ...> of field variables.
 *  Only field variables which are defined on the specified mesh will be output.
 *  The BasisOnMeshType has to be the type of the field variables given in meshName.
 */
template<typename BasisOnMeshType, typename OutputFieldVariablesType>
class ExfileWriter
{
public:

  //! write exelem file to given stream, for the mesh with given meshName
  static void outputExelem(std::ostream &stream, OutputFieldVariablesType fieldVariables, std::string meshName, std::shared_ptr<BasisOnMeshType> mesh, int nFieldVariablesOfMesh);

  //! write exnode file to given stream, only output fieldVariables that are on a mesh with the given meshName 
  static void outputExnode(std::ostream &stream, OutputFieldVariablesType fieldVariables, std::string meshName, std::shared_ptr<BasisOnMeshType> mesh, int nFieldVariablesOfMesh);
};

// specialization for UnstructuredDeformable
template<int D, typename BasisFunctionType, typename OutputFieldVariablesType>
class ExfileWriter<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>, OutputFieldVariablesType>
{
public:
  typedef BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType> BasisOnMeshType;

  //! write exelem file to given stream, for the mesh with given meshName
  static void outputExelem(std::ostream &stream, OutputFieldVariablesType fieldVariables, std::string meshName, std::shared_ptr<BasisOnMeshType> mesh, int nFieldVariablesOfMesh);

  //! write exnode file to given stream, only output fieldVariables that are on a mesh with the given meshName 
  static void outputExnode(std::ostream &stream, OutputFieldVariablesType fieldVariables, std::string meshName, std::shared_ptr<BasisOnMeshType> mesh, int nFieldVariablesOfMesh);
};


};  // namespace

#include "output_writer/exfile/exfile_writer_structured.tpp"
#include "output_writer/exfile/exfile_writer_unstructured_deformable.tpp"
