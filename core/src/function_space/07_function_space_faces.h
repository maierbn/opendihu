#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "mesh/face_t.h"
#include "mesh/type_traits.h"
#include "function_space/06_function_space_dofs_nodes.h"

namespace FunctionSpace
{

/** Adds functionality to get dof numbers corresponding to the faces of the element.
 */
template<typename MeshType,typename BasisFunctionType, typename = MeshType>
class FunctionSpaceFaces
{};

/** Partial specialization for 1D meshes. In this case a face is a single point (1 dof)
 */
template<typename MeshType,typename BasisFunctionType>
class FunctionSpaceFaces<MeshType, BasisFunctionType, Mesh::isDim<1,MeshType>> :
  public FunctionSpaceDofsNodes<MeshType,BasisFunctionType>
{
public:

  //! inherit constructor
  using FunctionSpaceDofsNodes<MeshType,BasisFunctionType>::FunctionSpaceDofsNodes;

  //! get all dof indices of a face, note: dimension in FunctionSpaceBaseDim is current-1 (=0), in this case the dofIndices array has exactly 1 entry
  static void getFaceDofs(Mesh::face_t face, std::array<dof_no_t,FunctionSpaceBaseDim<0,BasisFunctionType>::nDofsPerElement()> &dofIndices);
};

/** Partial specialization for 2D meshes. In this case a face is a line
 */
template<typename MeshType,typename BasisFunctionType>
class FunctionSpaceFaces<MeshType, BasisFunctionType, Mesh::isDim<2,MeshType>> :
  public FunctionSpaceDofsNodes<MeshType,BasisFunctionType>
{
public:

  //! inherit constructor
  using FunctionSpaceDofsNodes<MeshType,BasisFunctionType>::FunctionSpaceDofsNodes;

  //! get all dof indices of a face, note: dimension in FunctionSpaceBaseDim is current-1 (=1)
  static void getFaceDofs(Mesh::face_t face, std::array<dof_no_t,FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerElement()> &dofIndices);
};

/** Partial specialization for 3D meshes. In this case a face is a real 2D face.
 */
template<typename MeshType,typename BasisFunctionType>
class FunctionSpaceFaces<MeshType, BasisFunctionType, Mesh::isDim<3,MeshType>> :
  public FunctionSpaceDofsNodes<MeshType,BasisFunctionType>
{
public:

  //! inherit constructor
  using FunctionSpaceDofsNodes<MeshType,BasisFunctionType>::FunctionSpaceDofsNodes;

  //! get all dof indices of a face, note: dimension in FunctionSpaceBaseDim is current-1 (=2)
  static void getFaceDofs(Mesh::face_t face, std::array<dof_no_t,FunctionSpaceBaseDim<2,BasisFunctionType>::nDofsPerElement()> &dofIndices);
};

}  // namespace

#include "function_space/07_function_space_faces.tpp"
