#pragma once

#include <Python.h>  // has to be the first included header
#include <iostream>
#include <vector>
#include <memory>

#include "basis_on_mesh/05_basis_on_mesh.h"

namespace OutputWriter
{   

/** base class of ExfileWriter that writes exelem and exnode files of given field variables
 */
template<typename BasisOnMeshType>
class ExfileWriter
{};


// specialization for RegularFixed
template<int D, typename BasisFunctionType>
class ExfileWriter<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>
{
public:
  typedef BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType> BasisOnMeshType;
 
  //! write exnode file to given stream
  static void outputExelem(std::ostream &stream, std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> fieldVariables);
  
  //! write exnode file to given stream
  static void outputExnode(std::ostream &stream, std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> fieldVariables);
};


// specialization for StructuredDeformable 
template<int D, typename BasisFunctionType>
class ExfileWriter<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>
{
public:
  typedef BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType> BasisOnMeshType;
  
  //! write exnode file to given stream
  static void outputExelem(std::ostream &stream, std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> fieldVariables);
  
  //! write exnode file to given stream
  static void outputExnode(std::ostream &stream, std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> fieldVariables);
};


// specialization for UnstructuredDeformable
template<int D, typename BasisFunctionType>
class ExfileWriter<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>
{
public:
  typedef BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType> BasisOnMeshType;
  
  //! write exnode file to given stream
  static void outputExelem(std::ostream &stream, std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> fieldVariables);
  
  //! write exnode file to given stream
  static void outputExnode(std::ostream &stream, std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> fieldVariables);
};

};  // namespace

#include "output_writer/exfile/exfile_writer_regular_fixed.tpp"
#include "output_writer/exfile/exfile_writer_structured_deformable.tpp"
#include "output_writer/exfile/exfile_writer_unstructured_deformable.tpp"
