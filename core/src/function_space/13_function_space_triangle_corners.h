#pragma once

#include <Python.h>  // has to be the first included header

#include "control/types.h"
#include "mesh/type_traits.h"
#include "function_space/12_function_space_xi.h"

namespace FunctionSpace
{

/** This class adds tetahedral elements in the corners to avoid near-singular elements in the muscle geometry
 */
template<typename MeshType,typename BasisFunctionType,typename Dummy=void>
class FunctionSpaceTriangleCorners :
  public FunctionSpacePointInElement<MeshType,BasisFunctionType>
{
public:

  //! inherit constructor
  using FunctionSpacePointInElement<MeshType,BasisFunctionType>::FunctionSpacePointInElement;

  //! evaluate the basis function corresponding to element-local dof dofIndex at xi, xi lives in [0,1]^D
  virtual double phi(int dofIndex, std::array<double,MeshType::dim()> xi, element_no_t elementNoLocal) const override;

  //! evaluate the derivative of Phi(xi) w.r.t xi_i, where i is given by derivativeIdx, i.e. Phi_{dofIndex,derivativeIdx}(xi)
  virtual double dphi_dxi(int dofIndex, int derivativeIdx, std::array<double,MeshType::dim()> xi, element_no_t elementNoLocal) const override;

  //! set the dependent dofs in the given field variable by interpolating the independent dofs of the triangle basis, this implementation does nothing
  void interpolateNonDofValuesInFieldVariable(
    std::shared_ptr<FieldVariable::FieldVariableBaseFunctionSpace<FunctionSpace<MeshType,BasisFunctionType>>> fieldVariable,
    int componentNo) const{}

  //! set the dependent dofs in the given field variable by interpolating the independent dofs of the triangle basis
  template<int nComponents>
  void interpolateNonDofValuesInFieldVariable(
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace<::Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<2>>, nComponents>> fieldVariable) const{}
};

/** This class adds tetahedral elements in the corners to avoid near-singular elements in the muscle geometry
 * Partial specialization for 3D linear structured deformable mesh
 */
template<typename Dummy>
class FunctionSpaceTriangleCorners<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<1>,Dummy> :
  public FunctionSpacePointInElement<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<1>>
{
public:

  //! inherit constructor
  using FunctionSpacePointInElement<::Mesh::StructuredDeformableOfDimension<3>, ::BasisFunction::LagrangeOfOrder<1>>::FunctionSpacePointInElement;

  //! evaluate the basis function corresponding to element-local dof dofIndex at xi, xi lives in [0,1]^D
  double phi(int dofIndex, std::array<double,3> xi, element_no_t elementNoLocal) const;

  //! evaluate the derivative of Phi(xi) w.r.t xi_i, where i is given by derivativeIdx, i.e. Phi_{dofIndex,derivativeIdx}(xi)
  double dphi_dxi(int dofIndex, int derivativeIdx, std::array<double,3> xi, element_no_t elementNoLocal) const;

  //! set the dependent dofs in the given field variable by interpolating the independent dofs of the triangle basis
  void interpolateNonDofValuesInFieldVariable(
    std::shared_ptr<FieldVariable::FieldVariableBaseFunctionSpace<FunctionSpace<::Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<1>>>> fieldVariable,
    int componentNo) const;

  //! set the dependent dofs in the given field variable by interpolating the independent dofs of the triangle basis
  template<int nComponents>
  void interpolateNonDofValuesInFieldVariable(
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace<::Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<1>>,nComponents>> fieldVariable) const;
};

/** This class adds tetahedral elements in the corners to avoid near-singular elements in the muscle geometry
 * Partial specialization for 3D quadratic structured deformable mesh
 */
template<typename Dummy>
class FunctionSpaceTriangleCorners<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<2>, Dummy> :
  public FunctionSpacePointInElement<Mesh::StructuredDeformableOfDimension<3>, BasisFunction::LagrangeOfOrder<2>>
{
public:

  //! inherit constructor
  using FunctionSpacePointInElement<::Mesh::StructuredDeformableOfDimension<3>, ::BasisFunction::LagrangeOfOrder<2>>::FunctionSpacePointInElement;

  //! evaluate the basis function corresponding to element-local dof dofIndex at xi, xi lives in [0,1]^D
  double phi(int dofIndex, std::array<double,3> xi, element_no_t elementNoLocal) const;

  //! evaluate the derivative of Phi(xi) w.r.t xi_i, where i is given by derivativeIdx, i.e. Phi_{dofIndex,derivativeIdx}(xi)
  double dphi_dxi(int dofIndex, int derivativeIdx, std::array<double,3> xi, element_no_t elementNoLocal) const;

  //! set the dependent dofs in the given field variable by interpolating the independent dofs of the triangle basis
  void interpolateNonDofValuesInFieldVariable(
    std::shared_ptr<FieldVariable::FieldVariableBaseFunctionSpace<FunctionSpace<::Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<2>>>> fieldVariable,
    int componentNo) const;

  //! set the dependent dofs in the given field variable by interpolating the independent dofs of the triangle basis
  template<int nComponents>
  void interpolateNonDofValuesInFieldVariable(
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace<::Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<2>>,nComponents>> fieldVariable) const;
};

}  // namespace

#include "function_space/13_function_space_triangle_corners.tpp"
#include "function_space/13_function_space_triangle_corners_linear.tpp"
#include "function_space/13_function_space_triangle_corners_quadratic.tpp"
