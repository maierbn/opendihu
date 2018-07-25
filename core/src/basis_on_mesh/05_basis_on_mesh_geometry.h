#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"

#include "basis_on_mesh/04_basis_on_mesh_numbers_structured.h"
#include "basis_on_mesh/04_basis_on_mesh_data_unstructured.h"

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

  //! return a field variable with given name, this is not implemented for structured meshes since there are no extra stored field variables, only for unstructured meshes is it implemented and then stores field variables that were present in parsed exfiles.
  std::shared_ptr<FieldVariableBaseType> fieldVariable(std::string name);

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

protected:
 
  //! set up the geometry field
  virtual void initializeGeometryField() = 0;
  
  //! set the values of the geometry field
  virtual void setGeometryFieldValues() = 0;
  
};

}  // namespace

#include "basis_on_mesh/05_basis_on_mesh_geometry.tpp"
#include "basis_on_mesh/05_basis_on_mesh_geometry_structured.tpp"
#include "basis_on_mesh/05_basis_on_mesh_geometry_unstructured.tpp"
