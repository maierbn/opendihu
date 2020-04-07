#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"

#include "function_space/01_function_space_function.h"
#include "basis_function/lagrange.h"
#include "mesh/type_traits.h"

namespace FunctionSpace
{

/** Class with general algorithm to compute jacobian from basis functions.
 *  Note that this makes no sense for complete polynomials because these are not used to describe geometry.
 */
template<typename MeshType,typename BasisFunctionType,typename dummy = MeshType>
class FunctionSpaceJacobian :
  public FunctionSpaceFunction<MeshType,BasisFunctionType>
{
public:
  //! inherit constructor
  using FunctionSpaceFunction<MeshType,BasisFunctionType>::FunctionSpaceFunction;

  //! compute the (geometry) jacobian matrix, geometryField is the node positions for Lagrange basis, node positions and derivatives for Hermite basis
  template<typename Vec3>
  static std::array<Vec3,MeshType::dim()> computeJacobian(const std::array<Vec3,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> &geometryField,
                                                          const std::array<double,MeshType::dim()> xi);
};

/** partial specialization for linear Lagrange, D=1
 */
template<typename MeshType>
class FunctionSpaceJacobian<MeshType,BasisFunction::LagrangeOfOrder<1>,Mesh::isDim<1,MeshType>> :
  public FunctionSpaceFunction<MeshType,BasisFunction::LagrangeOfOrder<1>>
{
public:
  //! inherit constructor
  using FunctionSpaceFunction<MeshType,BasisFunction::LagrangeOfOrder<1>>::FunctionSpaceFunction;

  //! compute the jacobian matrix, geometryField is the node positions for Lagrange basis, node positions and derivatives for Hermite basis
  template<typename Vec3>
  static std::array<Vec3,1> computeJacobian(const std::array<Vec3,FunctionSpaceFunction<MeshType,BasisFunction::LagrangeOfOrder<1>>::nDofsPerElement()> &geometryField,
                                            const std::array<double,1> xi);

};

// partial specialization for linear Lagrange, D=2
template<typename MeshType>
class FunctionSpaceJacobian<MeshType,BasisFunction::LagrangeOfOrder<1>,Mesh::isDim<2,MeshType>> :
  public FunctionSpaceFunction<MeshType,BasisFunction::LagrangeOfOrder<1>>
{
public:
  //! inherit constructor
  using FunctionSpaceFunction<MeshType,BasisFunction::LagrangeOfOrder<1>>::FunctionSpaceFunction;

  //! compute the jacobian matrix, geometryField is the node positions for Lagrange basis, node positions and derivatives for Hermite basis
  template<typename Vec3>
  static std::array<Vec3,2> computeJacobian(const std::array<Vec3,FunctionSpaceFunction<MeshType,BasisFunction::LagrangeOfOrder<1>>::nDofsPerElement()> &geometryField,
                                            const std::array<double,2> xi);

};

// partial specialization for linear Lagrange, D=3
template<typename MeshType>
class FunctionSpaceJacobian<MeshType,BasisFunction::LagrangeOfOrder<1>,Mesh::isDim<3,MeshType>> :
  public FunctionSpaceFunction<MeshType,BasisFunction::LagrangeOfOrder<1>>
{
public:
  //! inherit constructor
  using FunctionSpaceFunction<MeshType,BasisFunction::LagrangeOfOrder<1>>::FunctionSpaceFunction;

  //! compute the jacobian matrix, geometryField is the node positions for Lagrange basis, node positions and derivatives for Hermite basis
  template<typename Vec3>
  static std::array<Vec3,3> computeJacobian(const std::array<Vec3,FunctionSpaceFunction<MeshType,BasisFunction::LagrangeOfOrder<1>>::nDofsPerElement()> &geometryField,
                                            const std::array<double,3> xi);

};

}  // namespace

#include "function_space/02_function_space_jacobian.tpp"
