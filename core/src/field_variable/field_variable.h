#pragma once

#include "field_variable/structured/06_field_variable_set_get_structured_deformable.h"
#include "field_variable/structured/07_field_variable_set_get_component_dependent_structured_regular_fixed.h"
#include "field_variable/unstructured/04_field_variable_set_get_component_dependent_unstructured_deformable.h"

namespace FieldVariable
{

template<typename BasisOnMeshType,int nComponents>
class FieldVariable :
  public FieldVariableSetGet<BasisOnMeshType,nComponents>
{ 
public:
  //! inherited constructor 
//   typedef BasisOnMesh::BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType> BasisOnMeshType;
  using FieldVariableSetGet<BasisOnMeshType,nComponents>::FieldVariableSetGet;
  
  typedef BasisOnMeshType BasisOnMesh;
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
template<typename BasisOnMeshType,int nComponents>
std::ostream &operator<<(std::ostream &stream, const FieldVariable<BasisOnMeshType,nComponents> &rhs)
{
  rhs.output(stream);
  return stream;
}

};  // namespace
