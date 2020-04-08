#pragma once

#include <Python.h>  // has to be the first included header

#include "function_space/08_function_space_nodes.h"

namespace FunctionSpace
{

//! forward declaration for field variables
template<typename MeshType,typename BasisFunctionType>
class FunctionSpace;


/** Class for structured meshes
 */
template<typename MeshType,typename BasisFunctionType>
class FunctionSpaceStructuredFindPositionBase :
  public FunctionSpaceNodes<MeshType,BasisFunctionType>
{
public:

  //! inherit constructor
  using FunctionSpaceNodes<MeshType,BasisFunctionType>::FunctionSpaceNodes;

  //! store a ghost mesh which is a neighouring mesh with only one layer of elements, this will be used by pointIsInElement and findPosition
  void setGhostMesh(Mesh::face_t face, const std::shared_ptr<FunctionSpace<MeshType,BasisFunctionType>> ghostMesh);

  //! get the element no and the xi value of the point, return true if the point is inside the mesh or false otherwise. Start search at given elementNo
  //! ghostMeshNo: -1 means main mesh, 0-5 means ghost Mesh with respecitve Mesh::face_t
  bool findPosition(Vec3 point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,MeshType::dim()> &xi, bool startSearchInCurrentElement, double &residual, double xiTolerance = 1e-4);

  //! check if the point lies inside the element, if yes, return true and set xi to the value of the point, defined in 11_function_space_xi.h
  virtual bool pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,MeshType::dim()> &xi, double &residual, double xiTolerance) = 0;

  //! print via VLOG(1) << which ghostMesh_ variables are set
  void debugOutputGhostMeshSet();

protected:

  //! check if the point is in a neighbouring element to elementNo on ghostMeshNo (-1=main mesh, 0-5=ghost mesh on respective face, 0=face0Minus, 1=face0Plus, etc.), return true if the element was found amoung the neighbours
  //! set elementNo, ghostMeshNo and xi appropriately
  virtual bool checkNeighbouringElements(const Vec3 &point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,MeshType::dim()> &xi, double &residual) = 0;

  std::array<std::shared_ptr<FunctionSpace<MeshType,BasisFunctionType>>,6> ghostMesh_;   // neighbouring functionSpaces of the local domain, i.e. containing ghost elements, this is used by findPosition
};

}  // namespace

#include "function_space/09_function_space_find_position_base.tpp"
