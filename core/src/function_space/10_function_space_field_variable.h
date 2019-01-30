#pragma once

#include <Python.h>  // has to be the first included header

#include "function_space/09_function_space_find_position.h"

namespace FunctionSpace
{

//! forward declaration for field variables
template<typename MeshType,typename BasisFunctionType>
class FunctionSpace;

/** FunctionSpace derives from the mesh and adds functionality to get dof and node numbers, phi and gradient.
 */
template<typename MeshType,typename BasisFunctionType>
class FunctionSpaceFieldVariable :
  public FunctionSpaceFindPosition<MeshType,BasisFunctionType>
{
public:

  //! inherit constructor
  using FunctionSpaceFindPosition<MeshType,BasisFunctionType>::FunctionSpaceFindPosition;

  //! create a non-geometry field field variable with no values being set, with given component names
  std::shared_ptr<FieldVariable::FieldVariableBaseFunctionSpace<FunctionSpace<MeshType,BasisFunctionType>>> createFieldVariable(std::string name, std::vector<std::string> componentNames);

  //! create a non-geometry field field variable with no values being set, with given number of components, the component names will be the numbers
  std::shared_ptr<FieldVariable::FieldVariableBaseFunctionSpace<FunctionSpace<MeshType,BasisFunctionType>>> createFieldVariable(std::string name, int nComponents=1);

  //! create a non-geometry field field variable with no values being set, with given number of components, the component names will be the numbers
  template <int nComponents>
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace<MeshType,BasisFunctionType>,nComponents>> createFieldVariable(std::string name);

  //! create a non-geometry field field variable with no values being set, with given number of components and component names
  template <int nComponents>
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace<MeshType,BasisFunctionType>,nComponents>> createFieldVariable(std::string name, std::vector<std::string> componentNames);

  //! get the number of scale factors that are stored for the global element no
  int getNumberScaleFactors(element_no_t elementGlobalNo);
  
  //! return an array of the gradients of all nodal basis functions, evaluated at xi
  std::array<std::array<double,MeshType::dim()>,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()>
  getGradPhi(std::array<double,MeshType::dim()> xi) const;

  //! interpolate the nComponents values within an element at the given xi position using the basis functions
  template <int nComponents>
  std::array<double,nComponents> interpolateValueInElement(std::array<std::array<double,nComponents>,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> &elementalDofValues,
                                                  std::array<double,MeshType::dim()> xi) const;

  //! interpolate the value within an element at the given xi position using the basis functions
  double interpolateValueInElement(std::array<double,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> &elementalDofValues,
                                   std::array<double,MeshType::dim()> xi) const;

  //! interpolate the gradient of a scalar field within an element at the given xi position using the basis functions
  //! the inverseJacobianParameterSpace can be computed by getInverseJacobian
  std::array<double,MeshType::dim()> interpolateGradientInElement(std::array<double,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> &elementalDofValues,
                                                                  Tensor2<MeshType::dim()> inverseJacobianParameterSpace, std::array<double,MeshType::dim()> xi) const;

  //! compute the normal in world space, normal to face at xi, use the given geometry values, that can by obtained by fieldVariable->getElementValues(elementNo, geometryValues) or mesh->getElementGeometry(elementNo, geometryValues)
  Vec3 getNormal(Mesh::face_t face, std::array<Vec3,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> geometryValues, std::array<double,MeshType::dim()> xi);

  //! compute the normal in world space, normal to face at xi
  Vec3 getNormal(Mesh::face_t face, element_no_t elementNo, std::array<double,MeshType::dim()> xi);

  //! Compute the inverseJacobian that is needed to transform a gradient vector from parameter space to world space, for an element at a xi position.
  //! This version of the method needs the values of the geometry field, if the jacobian is needed at multiple positions in the same element, these values can be retrieved once and used for all computations of the jacobians.
  //! There is also the convienience method which does not need the geometryValues but gets them itself.
  //! The following properties of the jacobian hold:
  //! jacobianParameterSpace[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx
  //! inverseJacobianParameterSpace[columnIdx][rowIdx] = dxi_rowIdx/dX_columnIdx because of inverse function theorem
  Tensor2<MeshType::dim()> getInverseJacobian(std::array<Vec3,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> &geometryValues, element_no_t elementNo, std::array<double,MeshType::dim()> xi);
  
  //! Compute the inverseJacobian that is needed to transform a gradient vector from parameter space to world space, for an element at a xi position.
  //! The following properties of the jacobian hold:
  //! jacobianParameterSpace[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx
  //! inverseJacobianParameterSpace[columnIdx][rowIdx] = dxi_rowIdx/dX_columnIdx because of inverse function theorem
  Tensor2<MeshType::dim()> getInverseJacobian(element_no_t elementNo, std::array<double,MeshType::dim()> xi);
};

}  // namespace

#include "function_space/10_function_space_field_variable.tpp"
