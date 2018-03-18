#pragma once

#include "field_variable/structured/04_field_variable_set_get_structured_deformable.h"
#include "field_variable/structured/04_field_variable_set_get_structured_regular_fixed.h"
#include "field_variable/unstructured/02_field_variable_set_get_unstructured_deformable.h"

namespace FieldVariable
{

template<typename BasisOnMeshType>
class FieldVariable :
  public FieldVariableSetGet<BasisOnMeshType>
{ 
public:
  //! inherited constructor 
  using FieldVariableSetGet<BasisOnMeshType>::FieldVariableSetGet;
};
// /** FieldVariable class for StructuredRegularFixedOfDimension mesh
//  */
// template<int D, typename BasisFunctionType>
// class FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>> :
//   public FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>,
//   public Interface<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>
// { 
// public:
//   //! inherited constructor 
//   using FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>::FieldVariableSetGet;
//   
//   typedef BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType> BasisOnMeshType;
// };
// 
// /** FieldVariable class for StructuredDeformableOfDimension mesh
//  */
// template<int D, typename BasisFunctionType>
// class FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>> :
//   public FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>,
//   public Interface<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>
// { 
// public:
//   //! inherited constructor 
//   using FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>::FieldVariableSetGet;
//   
//   typedef BasisOnMesh::BasisOnMesh<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType> BasisOnMeshType;
// };
// 
// /** FieldVariable class for UnstructuredDeformable mesh
//  */
// template<int D, typename BasisFunctionType>
// class FieldVariable<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>> :
//   public FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>,
//   public Interface<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>
// { 
// public:
//   //! inherited constructor 
//   using FieldVariableSetGet<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>::FieldVariableSetGet;
//   
//   typedef BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType> BasisOnMeshType;
// };

// output operator
template<typename BasisOnMeshType>
std::ostream &operator<<(std::ostream &stream, const FieldVariable<BasisOnMeshType> &rhs)
{
  rhs.output(stream);
  return stream;
}

};  // namespace
