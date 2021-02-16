#pragma once

#include <Python.h>  // has to be the first included header

#include "function_space/08_function_space_nodes.h"
#include "function_space/09_function_space_structured_check_neighbouring_elements.h"

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
  //! inherit constructor
  using FunctionSpaceNodes<MeshType,BasisFunctionType>::FunctionSpaceNodes;

};

/** Partial specialization for structured meshes
 */
template<typename MeshType,typename BasisFunctionType>
class FunctionSpaceFindPosition<MeshType,BasisFunctionType,Mesh::isStructured<MeshType>> :
  public FunctionSpaceStructuredCheckNeighbouringElements<MeshType,BasisFunctionType>
{
public:

  //! inherit constructor
  using FunctionSpaceStructuredCheckNeighbouringElements<MeshType,BasisFunctionType>::FunctionSpaceStructuredCheckNeighbouringElements;
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
  bool findPosition(Vec3 point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,D> &xi, bool startSearchInCurrentElement, double &residual, bool &searchedAllElements, double xiTolerance = 1e-4);

  //! check if the point lies inside the element, if yes, return true and set xi to the value of the point, defined in 11_function_space_xi.h
  virtual bool pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,D> &xi, double &residual, double xiTolerance) = 0;

  //! return a nullptr,  for structured meshes this is a pointer to the ghost mesh indexed by faceOrEdge
  std::shared_ptr<FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>> ghostMesh(Mesh::face_or_edge_t faceOrEdge);
};

/** Partial specialization for composite meshes
 */
template<int D,typename BasisFunctionType>
class FunctionSpaceFindPosition<Mesh::CompositeOfDimension<D>,BasisFunctionType,Mesh::CompositeOfDimension<D>> :
  public FunctionSpaceStructuredCheckNeighbouringElements<Mesh::CompositeOfDimension<D>,BasisFunctionType>
{
public:

  //! inherit constructor
  using FunctionSpaceStructuredCheckNeighbouringElements<Mesh::CompositeOfDimension<D>,BasisFunctionType>::FunctionSpaceStructuredCheckNeighbouringElements;
};


}  // namespace

#include "function_space/09_function_space_find_position.tpp"
