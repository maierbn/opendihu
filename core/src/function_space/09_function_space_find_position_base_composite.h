#pragma once

#include <Python.h>  // has to be the first included header

#include "function_space/08_function_space_nodes.h"

namespace FunctionSpace
{

//! forward declaration for field variables
template<typename MeshType,typename BasisFunctionType>
class FunctionSpace;

/** partial specialization for composite mesh
 */
template<int D,typename BasisFunctionType>
class FunctionSpaceStructuredFindPositionBase<Mesh::CompositeOfDimension<D>,BasisFunctionType> :
  public FunctionSpaceNodes<Mesh::CompositeOfDimension<D>,BasisFunctionType>
{
public:

  typedef Mesh::CompositeOfDimension<D> MeshType;

  //! inherit constructor
  using FunctionSpaceNodes<MeshType,BasisFunctionType>::FunctionSpaceNodes;

  //! store a ghost mesh which is a neighouring mesh with only one layer of elements, this will be used by pointIsInElement and findPosition
  void setGhostMesh(Mesh::face_t face, const std::shared_ptr<FunctionSpace<MeshType,BasisFunctionType>> ghostMesh){}

  //! get the element no and the xi value of the point, return true if the point is inside the mesh or false otherwise. Start search at given elementNo
  //! ghostMeshNo: -1 means main mesh, 0-5 means ghost Mesh with respective Mesh::face_t
  bool findPosition(Vec3 point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,MeshType::dim()> &xi, bool startSearchInCurrentElement, double &residual, bool &searchedAllElements, double xiTolerance = 1e-4);

  //! check if the point lies inside the element, if yes, return true and set xi to the value of the point, defined in 11_function_space_xi.h
  virtual bool pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,MeshType::dim()> &xi, double &residual, double xiTolerance) = 0;

  //! print via VLOG(1) << which ghostMesh_ variables are set
  void debugOutputGhostMeshSet(){}

  //! return the sub mesh no. where the last point was found by findPosition
  int subMeshNoWherePointWasFound();

protected:
  int subMeshNoWherePointWasFound_ = 0;    //< findPositions sets this to the subMeshNo in which the point was found
};

}  // namespace

#include "function_space/09_function_space_find_position_base_composite.tpp"
