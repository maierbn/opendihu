#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "mesh/face_t.h"
#include "mesh/type_traits.h"
#include "function_space/10_function_space_field_variable.h"

namespace FunctionSpace
{

// forward declarations (declarations at the end of this file)
template<typename MeshType,typename BasisFunctionType, typename DummyForTraits = MeshType>
class ComputeXiApproximation;

template<typename MeshType>
class ComputeXiApproximation<MeshType, BasisFunction::LagrangeOfOrder<1>, Mesh::isDeformableWithDim<3,MeshType>>;

/** Helper class shared by all specializations that checks if point lies inside the bounding box of an element
 */
template<typename MeshType,typename BasisFunctionType>
class FunctionSpaceXi :
  public FunctionSpaceFieldVariable<MeshType,BasisFunctionType>
{
public:
  //! inherit constructor
  using FunctionSpaceFieldVariable<MeshType,BasisFunctionType>::FunctionSpaceFieldVariable;

  //! check if the point lies outside the element with given geometryValues, if yes, return true, if it returns false this does not mean the point has to lie inside
  bool pointIsOutsideBoundingBox(Vec3 point, const std::array<Vec3,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> &geometryValues, double xiTolerance) const;
  
  //! check if the point is one of the node positions given in geometryValues, if so, return true and set xi accordingly
  bool pointIsNodePosition(Vec3 point, const std::array<Vec3,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> &geometryValues, std::array<double,MeshType::dim()> &xi) const;
};

// -----------------------------------------------------

/** Adds functionality to get the xi value, i.e. the local coordinate, for a point inside an element
 */
template<typename MeshType,typename BasisFunctionType, typename = MeshType>
class FunctionSpacePointInElement :
  public ComputeXiApproximation<MeshType,BasisFunctionType>
{
public:
  //! inherit constructor
  using ComputeXiApproximation<MeshType,BasisFunctionType>::ComputeXiApproximation;

  //! check if the point lies inside the element, if yes, return true and set xi to the value of the point
  bool pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,MeshType::dim()> &xi, double &residual, double xiTolerance = 1e-4);
};

/** Partial specialization for 1D StructuredRegularFixed meshes
 */
template<typename BasisFunctionType>
class FunctionSpacePointInElement<Mesh::StructuredRegularFixedOfDimension<1>, BasisFunctionType> :
  public FunctionSpaceFieldVariable<Mesh::StructuredRegularFixedOfDimension<1>,BasisFunctionType>
{
public:

  //! inherit constructor
  using FunctionSpaceFieldVariable<Mesh::StructuredRegularFixedOfDimension<1>,BasisFunctionType>::FunctionSpaceFieldVariable;

  //! check if the point lies inside the element, if yes, return true and set xi to the value of the point
  bool pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,1> &xi, double &residual, double xiTolerance = 1e-4);
};

/** Partial specialization for 2D StructuredRegularFixed meshes
 */
template<typename BasisFunctionType>
class FunctionSpacePointInElement<Mesh::StructuredRegularFixedOfDimension<2>, BasisFunctionType> :
  public FunctionSpaceFieldVariable<Mesh::StructuredRegularFixedOfDimension<2>,BasisFunctionType>
{
public:

  //! inherit constructor
  using FunctionSpaceFieldVariable<Mesh::StructuredRegularFixedOfDimension<2>,BasisFunctionType>::FunctionSpaceFieldVariable;

  //! check if the point lies inside the element, if yes, return true and set xi to the value of the point
  bool pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,2> &xi, double &residual, double xiTolerance = 1e-4);
};

/** Partial specialization for 3D StructuredRegularFixed meshes
 */
template<typename BasisFunctionType>
class FunctionSpacePointInElement<Mesh::StructuredRegularFixedOfDimension<3>, BasisFunctionType> :
  public FunctionSpaceFieldVariable<Mesh::StructuredRegularFixedOfDimension<3>,BasisFunctionType>
{
public:

  //! inherit constructor
  using FunctionSpaceFieldVariable<Mesh::StructuredRegularFixedOfDimension<3>,BasisFunctionType>::FunctionSpaceFieldVariable;

  //! check if the point lies inside the element, if yes, return true and set xi to the value of the point
  bool pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,3> &xi, double &residual, double xiTolerance = 1e-4);
};

/** Partial specialization for 1D deformable meshes and linear shape functions
 */
template<typename MeshType>
class FunctionSpacePointInElement<MeshType, BasisFunction::LagrangeOfOrder<1>, Mesh::isDeformableWithDim<1,MeshType>> :
  public FunctionSpaceFieldVariable<MeshType,BasisFunction::LagrangeOfOrder<1>>
{
public:

  //! inherit constructor
  using FunctionSpaceFieldVariable<MeshType,BasisFunction::LagrangeOfOrder<1>>::FunctionSpaceFieldVariable;

  //! check if the point lies inside the element, if yes, return true and set xi to the value of the point
  bool pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,1> &xi, double &residual, double xiTolerance = 1e-4);
};

/** Partial specialization for 2D deformable meshes and linear shape functions
 */

template<typename MeshType>
class FunctionSpacePointInElement<MeshType, BasisFunction::LagrangeOfOrder<1>, Mesh::isDeformableWithDim<2,MeshType>> :
  public FunctionSpaceFieldVariable<MeshType,BasisFunction::LagrangeOfOrder<1>>
{
public:

  //! inherit constructor
  using FunctionSpaceFieldVariable<MeshType,BasisFunction::LagrangeOfOrder<1>>::FunctionSpaceFieldVariable;

  //! check if the point lies inside the element, if yes, return true and set xi to the value of the point
  bool pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,2> &xi, double &residual, double xiTolerance = 1e-4);
};

// --------------------------------------------------

/** Helper class that provides an initial guess for xi from point as initial value for the newton scheme
 * In this general class the initial guess is 0.
 */
template<typename MeshType,typename BasisFunctionType, typename DummyForTraits>
class ComputeXiApproximation :
  public FunctionSpaceXi<MeshType,BasisFunctionType>
{
public:
 
  //! inherit constructor
  using FunctionSpaceXi<MeshType,BasisFunctionType>::FunctionSpaceXi;

protected:
  //! compute an approximation for the xi values of point, assuming it is inside the element
  void computeApproximateXiForPoint(Vec3 point, element_no_t elementNo, std::array<double,MeshType::dim()> &xi);

};

/** Specialization for 3D and linear Lagrange geometry, compute initial xi by a fast heuristic.
 */
template<typename MeshType>
class ComputeXiApproximation<MeshType, BasisFunction::LagrangeOfOrder<1>, Mesh::isDeformableWithDim<3,MeshType>> :
  public FunctionSpaceXi<MeshType,BasisFunction::LagrangeOfOrder<1>>
{
public:
 
  //! inherit constructor
  using FunctionSpaceXi<MeshType,BasisFunction::LagrangeOfOrder<1>>::FunctionSpaceXi;

protected:
 
  //! check if the point lies inside the space which is spanned by a 4-point tetrahedron, computes the barycentric coordinates of the point 
  //! the return value is whether xi1,xi2,xi3 >= 0. correctOrientation specifies if  the components of xi need to be reversed as 1-xi
  bool pointIsInTetrahedron(Vec3 point, std::array<Vec3,4> tetrahedron, std::array<bool,3> correctOrientation, std::array<double,3> &xi);
  
  //! check if the point is inside the element via a fast heuristic
  bool pointIsInElementQuick(Vec3 point, const std::array<Vec3, FunctionSpaceBaseDim<3,BasisFunction::LagrangeOfOrder<1>>::nDofsPerElement()> &geometryValues) const;
   
  //! compute an approximation for the xi values of point, assuming it is inside the element
  void computeApproximateXiForPoint(Vec3 point, element_no_t elementNo, std::array<double,3> &xi);

}; 

}  // namespace

#include "function_space/11_function_space_xi.tpp"
