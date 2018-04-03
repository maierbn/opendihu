#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"

#include "basis_on_mesh/03_basis_on_mesh_numbers_structured.h"
#include "basis_on_mesh/03_basis_on_mesh_data_unstructured.h"

namespace BasisOnMesh
{

/** base class for structured meshes, geometry field is declared here structured meshes
 */
template<typename MeshType,typename BasisFunctionType>
class BasisOnMeshGeometryData :
  public BasisOnMeshNumbers<MeshType, BasisFunctionType>
{
public:
  //! inherit constructor
  using BasisOnMeshNumbers<MeshType,BasisFunctionType>::BasisOnMeshNumbers;
  
  typedef FieldVariable::FieldVariableBase<BasisOnMesh<MeshType,BasisFunctionType>> FieldVariableBaseType;  ///< the class typename of the a field variable
  typedef FieldVariable::FieldVariable<BasisOnMesh<MeshType,BasisFunctionType>,3> GeometryFieldType;  ///< the class typename of the geometry field variable

  //! return a field variable with given name, this is not implemented for structured meshes since there are no extra stored field variables, only for unstructured meshes is in implemented and then stores field variables that were present in parsed exfiles.
  std::shared_ptr<FieldVariableBaseType> fieldVariable(std::string name);
  /*{
    return nullptr;
  }*/

protected:
  std::unique_ptr<GeometryFieldType> geometryField_;     ///< the geometry field variable
  bool noGeometryField_ = false;                         ///< this is set if there is no geometry field stored. this is only needed for solid mechanics mixed formulation where the lower order basisOnMesh does not need its own geometry information
};

/** partial specialization for unstructured mesh, unstructured mesh already has geometry field declared inside BasisOnMeshDataUnstructured
 */
template<int D,typename BasisFunctionType>
class BasisOnMeshGeometryData<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType> :
  public BasisOnMeshDataUnstructured<D,BasisFunctionType>
{
public:
  //! inherited constructor 
  using BasisOnMeshDataUnstructured<D,BasisFunctionType>::BasisOnMeshDataUnstructured;
  
  typedef FieldVariable::FieldVariableBase<BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>> FieldVariableBaseType;  ///< the class typename of the a field variable
  
  //! return a field variable with given name, returns field variables that were present in parsed exfiles
  std::shared_ptr<FieldVariableBaseType> fieldVariable(std::string name)
  {
    if (this->fieldVariable_.find(name) != this->fieldVariable_.end())
      return this->fieldVariable_.at(name);
    else
      return nullptr;
  }

}; 

/** base class for all meshes, not complete polynomials as basis functions
 */
template<typename MeshType,typename BasisFunctionType,typename DummyForTraits=MeshType>
class BasisOnMeshGeometry :
  public BasisOnMeshGeometryData<MeshType,BasisFunctionType>,
  public std::enable_shared_from_this<BasisOnMeshGeometry<MeshType,BasisFunctionType>>
{
public:
  //! inherit constructor
  using BasisOnMeshGeometryData<MeshType,BasisFunctionType>::BasisOnMeshGeometryData;
  
  typedef FieldVariable::FieldVariable<BasisOnMesh<MeshType,BasisFunctionType>,3> GeometryFieldType;  ///< the class typename of the geometry field variable

  //! this assigns the geometry field variable's mesh pointer to this object, it is not possible from the constructor, therefore this extra method
  void initialize();
  
  //! return the geometry field entry (node position for Lagrange elements) of a specific dof
  Vec3 getGeometry(node_no_t dofGlobalNo) const;
  
  //! get all geometry entries for an element
  void getElementGeometry(element_no_t elementNo, std::array<Vec3, BasisOnMeshBaseDim<MeshType::dim(),BasisFunctionType>::nDofsPerElement()> &values);
  
  //! return the internal geometry field variable
  GeometryFieldType &geometryField();
  
  //! if the geometry field is set
  bool hasGeometryField();

};

/** partial specialization for structured mesh and complete polynomials
 */
/*
template<typename MeshType,int order>
class BasisOnMeshGeometry<MeshType,BasisFunction::CompletePolynomialOfDimensionAndOrder<MeshType::dim(),order>,Mesh::isDeformable<MeshType>> :
  public BasisOnMeshFunction<MeshType,BasisFunction::CompletePolynomialOfDimensionAndOrder<MeshType::dim(),order>>
{
public:
  //! inherit constructor
  using BasisOnMeshFunction<MeshType,BasisFunction::CompletePolynomialOfDimensionAndOrder<MeshType::dim(),order>>::BasisOnMeshJacobian;
  */
  /* TODO: remove
  //! return the global dof number of element-local dof dofIndex of element elementNo
  dof_no_t getDofNo(element_no_t elementNo, int dofIndex) const;
  
  //! return the global node number of element-local node nodeIndex of element elementNo
  node_no_t getNodeNo(element_no_t elementNo, int nodeIndex) const;
  
  //! get all dofs of a specific node
  void getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const;
  */
//};

}  // namespace

#include "basis_on_mesh/04_basis_on_mesh_geometry.tpp"
#include "basis_on_mesh/04_basis_on_mesh_geometry_structured.tpp"
#include "basis_on_mesh/04_basis_on_mesh_geometry_unstructured.tpp"
