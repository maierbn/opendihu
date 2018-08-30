#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"

#include "function_space/04_function_space_numbers_structured.h"
#include "function_space/04_function_space_data_unstructured.h"

namespace FunctionSpace
{

/** base class for structured meshes, geometry field is declared here structured meshes
 */
template<typename MeshType,typename BasisFunctionType>
class FunctionSpaceGeometryData :
  public FunctionSpaceNumbers<MeshType, BasisFunctionType>
{
public:
  //! inherit constructor
  using FunctionSpaceNumbers<MeshType,BasisFunctionType>::FunctionSpaceNumbers;

  typedef FieldVariable::FieldVariableBase<FunctionSpace<MeshType,BasisFunctionType>> FieldVariableBaseType;  ///< the class typename of the a field variable
  typedef FieldVariable::FieldVariable<FunctionSpace<MeshType,BasisFunctionType>,3> GeometryFieldType;  ///< the class typename of the geometry field variable

  //! return a field variable with given name, this is not implemented for structured meshes since there are no extra stored field variables, only for unstructured meshes is it implemented and then stores field variables that were present in parsed exfiles.
  std::shared_ptr<FieldVariableBaseType> fieldVariable(std::string name);

protected:
  std::unique_ptr<GeometryFieldType> geometryField_;     ///< the geometry field variable
  bool noGeometryField_ = false;                         ///< this is set if there is no geometry field stored. this is only needed for solid mechanics mixed formulation where the lower order basisOnMesh does not need its own geometry information
};

/** partial specialization for unstructured mesh, unstructured mesh already has geometry field declared inside FunctionSpaceDataUnstructured
 */
template<int D,typename BasisFunctionType>
class FunctionSpaceGeometryData<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType> :
  public FunctionSpaceDataUnstructured<D,BasisFunctionType>
{
public:
  //! inherited constructor
  using FunctionSpaceDataUnstructured<D,BasisFunctionType>::FunctionSpaceDataUnstructured;

  typedef FieldVariable::FieldVariableBase<FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>> FieldVariableBaseType;  ///< the class typename of the a field variable

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
class FunctionSpaceGeometry :
  public FunctionSpaceGeometryData<MeshType,BasisFunctionType>
{
public:
  //! inherit constructor
  using FunctionSpaceGeometryData<MeshType,BasisFunctionType>::FunctionSpaceGeometryData;

  typedef FieldVariable::FieldVariable<FunctionSpace<MeshType,BasisFunctionType>,3> GeometryFieldType;  ///< the class typename of the geometry field variable

  //! return the geometry field entry (node position for Lagrange elements) of a specific dof
  Vec3 getGeometry(node_no_t dofGlobalNo) const;

  //! get all geometry entries for an element
  void getElementGeometry(element_no_t elementNo, std::array<Vec3, FunctionSpaceBaseDim<MeshType::dim(),BasisFunctionType>::nDofsPerElement()> &values);

  //! return the internal geometry field variable
  GeometryFieldType &geometryField();

  //! if the geometry field is set
  bool hasGeometryField();

};

}  // namespace

#include "function_space/05_function_space_geometry.tpp"
#include "function_space/05_function_space_geometry_structured.tpp"
#include "function_space/05_function_space_geometry_unstructured.tpp"
