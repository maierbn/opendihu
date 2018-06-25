#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "mesh/face_t.h"
#include "mesh/type_traits.h"
#include "basis_on_mesh/09_basis_on_mesh_field_variable.h"

namespace BasisOnMesh
{

// forward declarations (declarations at the end of this file)
template<typename MeshType,typename BasisFunctionType, typename DummyForTraits = MeshType>
class ComputeXiApproximation;

template<typename MeshType>
class ComputeXiApproximation<MeshType, BasisFunction::LagrangeOfOrder<1>, Mesh::isDeformableWithDim<3,MeshType>>;

// -----------------------------------------------------

/** Adds functionality to get the xi value for a point inside an element
 */
template<typename MeshType,typename BasisFunctionType, typename = MeshType>
class BasisOnMeshXi : 
  public ComputeXiApproximation<MeshType,BasisFunctionType>
{
public:
  //! inherit constructor
  using ComputeXiApproximation<MeshType,BasisFunctionType>::ComputeXiApproximation;

  //! check if the point lies inside the element, if yes, return true and set xi to the value of the point
  bool pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,MeshType::dim()> &xi);
};

/** Partial specialization for 1D StructuredRegularFixed meshes
 */
template<typename BasisFunctionType>
class BasisOnMeshXi<Mesh::StructuredRegularFixedOfDimension<1>, BasisFunctionType> :
  public BasisOnMeshFieldVariable<Mesh::StructuredRegularFixedOfDimension<1>,BasisFunctionType>
{
public:

  //! inherit constructor
  using BasisOnMeshFieldVariable<Mesh::StructuredRegularFixedOfDimension<1>,BasisFunctionType>::BasisOnMeshFieldVariable;

  //! check if the point lies inside the element, if yes, return true and set xi to the value of the point
  bool pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,1> &xi);
};

/** Partial specialization for 2D StructuredRegularFixed meshes
 */
template<typename BasisFunctionType>
class BasisOnMeshXi<Mesh::StructuredRegularFixedOfDimension<2>, BasisFunctionType> :
  public BasisOnMeshFieldVariable<Mesh::StructuredRegularFixedOfDimension<2>,BasisFunctionType>
{
public:

  //! inherit constructor
  using BasisOnMeshFieldVariable<Mesh::StructuredRegularFixedOfDimension<2>,BasisFunctionType>::BasisOnMeshFieldVariable;

  //! check if the point lies inside the element, if yes, return true and set xi to the value of the point
  bool pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,2> &xi);
};

/** Partial specialization for 3D StructuredRegularFixed meshes
 */
template<typename BasisFunctionType>
class BasisOnMeshXi<Mesh::StructuredRegularFixedOfDimension<3>, BasisFunctionType> :
  public BasisOnMeshFieldVariable<Mesh::StructuredRegularFixedOfDimension<3>,BasisFunctionType>
{
public:

  //! inherit constructor
  using BasisOnMeshFieldVariable<Mesh::StructuredRegularFixedOfDimension<3>,BasisFunctionType>::BasisOnMeshFieldVariable;

  //! check if the point lies inside the element, if yes, return true and set xi to the value of the point
  bool pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,3> &xi);
};

/** Partial specialization for 1D deformable meshes and linear shape functions
 */
template<typename MeshType>
class BasisOnMeshXi<MeshType, BasisFunction::LagrangeOfOrder<1>, Mesh::isDeformableWithDim<1,MeshType>> :
  public BasisOnMeshFieldVariable<MeshType,BasisFunction::LagrangeOfOrder<1>>
{
public:

  //! inherit constructor
  using BasisOnMeshFieldVariable<MeshType,BasisFunction::LagrangeOfOrder<1>>::BasisOnMeshFieldVariable;

  //! check if the point lies inside the element, if yes, return true and set xi to the value of the point
  bool pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,1> &xi);
};

/** Partial specialization for 2D deformable meshes and linear shape functions
 */
template<typename MeshType>
class BasisOnMeshXi<MeshType, BasisFunction::LagrangeOfOrder<1>, Mesh::isDeformableWithDim<2,MeshType>> :
  public BasisOnMeshFieldVariable<MeshType,BasisFunction::LagrangeOfOrder<1>>
{
public:

  //! inherit constructor
  using BasisOnMeshFieldVariable<MeshType,BasisFunction::LagrangeOfOrder<1>>::BasisOnMeshFieldVariable;

  //! check if the point lies inside the element, if yes, return true and set xi to the value of the point
  bool pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,2> &xi);
};

// --------------------------------------------------

/** Helper class that provides an initial guess for xi from point as initial value for the newton scheme
 * In this general class the initial guess is 0.
 */
template<typename MeshType,typename BasisFunctionType, typename DummyForTraits>
class ComputeXiApproximation :
  public BasisOnMeshFieldVariable<MeshType,BasisFunctionType>
{
public:
 
  //! inherit constructor
  using BasisOnMeshFieldVariable<MeshType,BasisFunctionType>::BasisOnMeshFieldVariable;

protected:
  //! compute an approximation for the xi values of point, assuming it is inside the element
  void computeApproximateXiForPoint(Vec3 point, element_no_t elementNo, std::array<double,MeshType::dim()> &xi);

};

/** Specialization for 3D and linear Lagrange geometry, compute initial xi by a fast heuristic.
 */
template<typename MeshType>
class ComputeXiApproximation<MeshType, BasisFunction::LagrangeOfOrder<1>, Mesh::isDeformableWithDim<3,MeshType>> :
  public BasisOnMeshFieldVariable<MeshType,BasisFunction::LagrangeOfOrder<1>>
{
public:
 
  //! inherit constructor
  using BasisOnMeshFieldVariable<MeshType,BasisFunction::LagrangeOfOrder<1>>::BasisOnMeshFieldVariable;

protected:
 
  //! check if the point lies inside the space which is spanned by a 4-point tetrahedron, computes the barycentric coordinates of the point 
  //! the return value is whether xi1,xi2,xi3 >= 0. correctOrientation specifies if  the components of xi need to be reversed as 1-xi
  bool pointIsInTetrahedron(Vec3 point, std::array<Vec3,4> tetrahedron, std::array<bool,3> correctOrientation, std::array<double,3> &xi);
  
  //! check if the point is inside the element via a fast heuristic
  bool pointIsInElementQuick(Vec3 point, std::array<Vec3, BasisOnMeshBaseDim<3,BasisFunction::LagrangeOfOrder<1>>::nDofsPerElement()> &geometryValues);
   
  //! compute an approximation for the xi values of point, assuming it is inside the element
  void computeApproximateXiForPoint(Vec3 point, element_no_t elementNo, std::array<double,3> &xi);

}; 


 
}  // namespace

#include "basis_on_mesh/10_basis_on_mesh_xi.tpp"
