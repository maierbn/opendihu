#pragma once

#include <Python.h>  // has to be the first included header

#include "function_space/08_function_space_nodes.h"

namespace FunctionSpace
{

//! forward declaration for field variables
template<typename MeshType,typename BasisFunctionType>
class FunctionSpace;

/** FunctionSpace derives from the mesh and adds functionality to get dof and node numbers, phi and gradient.
 */
template<typename MeshType,typename BasisFunctionType,typename DummyForTraits=MeshType>
class FunctionSpaceFindPosition :
  public FunctionSpaceNodes<MeshType,BasisFunctionType>
{
};

/** Partial specialization for structured meshes
 */
template<typename MeshType,typename BasisFunctionType>
class FunctionSpaceFindPosition<MeshType,BasisFunctionType,Mesh::isStructured<MeshType>> :
  public FunctionSpaceNodes<MeshType,BasisFunctionType>
{
public:

  //! inherit constructor
  using FunctionSpaceNodes<MeshType,BasisFunctionType>::FunctionSpaceNodes;

  //! store a ghost mesh which is a neighouring mesh with only one layer of elements, this will be used by pointIsInElement and findPosition
  void setGhostMesh(Mesh::face_t face, const std::shared_ptr<FunctionSpace<MeshType,BasisFunctionType>> ghostMesh);

  //! get the element no and the xi value of the point, return true if the point is inside the mesh or false otherwise. Start search at given elementNo
  //! ghostMeshNo: -1 means main mesh, 0-5 means ghost Mesh with respecitve Mesh::face_t
  bool findPosition(Vec3 point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,MeshType::dim()> &xi, bool startSearchInCurrentElement);

  //! check if the point lies inside the element, if yes, return true and set xi to the value of the point, defined in 11_function_space_xi.h
  virtual bool pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,MeshType::dim()> &xi) = 0;

protected:

  //! check if the point is in a neighoubring element to elmentNo on ghostMeshNo (-1=main m0esh, 0-5=ghost mesh on respective face, 0=face0Minus, 1=face0Plus, etc.), return true if the element was found amoung the neighbours
  //! set elementNo, ghostMeshNo and xi appropriately
  bool checkNeighbouringElements(const Vec3 &point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,MeshType::dim()> &xi);

  int targetX_ = 0;   //! the x index of the neighbour (-1,0, or 1) where the last neighbouring element was found
  int targetY_ = 0;   //! the y index of the neighbour (-1,0, or 1) where the last neighbouring element was found
  int targetZ_ = 0;   //! the z index of the neighbour (-1,0, or 1) where the last neighbouring element was found

  std::array<std::shared_ptr<FunctionSpace<MeshType,BasisFunctionType>>,6> ghostMesh_;   // neighbouring functionSpaces of the local domain, i.e. containing ghost elements, this is used by findPosition
};

/** Partial specialization for unstructured meshes
 */
template<int D,typename BasisFunctionType>
class FunctionSpaceFindPosition<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType,Mesh::UnstructuredDeformableOfDimension<D>> :
  public FunctionSpaceNodes<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>
{
public:

  //! inherit constructor
  using FunctionSpaceNodes<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::FunctionSpaceNodes;

  //! get the element no and the xi value of the point, return true if the point is inside the mesh or false otherwise. Start search at given elementNo
  //! ghostMeshNo: -1 means main mesh, 0-5 means ghost Mesh with respecitve Mesh::face_t
  bool findPosition(Vec3 point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,D> &xi, bool startSearchInCurrentElement);

  //! check if the point lies inside the element, if yes, return true and set xi to the value of the point, defined in 11_function_space_xi.h
  virtual bool pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,D> &xi) = 0;

};



}  // namespace

#include "function_space/09_function_space_find_position.tpp"
