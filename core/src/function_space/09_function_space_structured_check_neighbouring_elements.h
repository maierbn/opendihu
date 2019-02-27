#pragma once

#include <Python.h>  // has to be the first included header

#include "function_space/09_function_space_find_position_base.h"

namespace FunctionSpace
{

//! forward declaration for field variables
template<typename MeshType,typename BasisFunctionType>
class FunctionSpace;

template<typename MeshType,typename BasisFunctionType, typename DummyForTraits=MeshType>
class FunctionSpaceStructuredCheckNeighbouringElements
{
};

/** Partial specialization for 1D meshes
 */
template<typename MeshType,typename BasisFunctionType>
class FunctionSpaceStructuredCheckNeighbouringElements<MeshType,BasisFunctionType,Mesh::isDim<1,MeshType>> :
  public FunctionSpaceStructuredFindPositionBase<MeshType,BasisFunctionType>
{
public:

  //! inherit constructor
  using FunctionSpaceStructuredFindPositionBase<MeshType,BasisFunctionType>::FunctionSpaceStructuredFindPositionBase;

protected:

  //! check if the point is in a neighbouring element to elementNo on ghostMeshNo (-1=main mesh, 0-5=ghost mesh on respective face, 0=face0Minus, 1=face0Plus, etc.), return true if the element was found amoung the neighbours
  //! set elementNo, ghostMeshNo and xi appropriately
  bool checkNeighbouringElements(const Vec3 &point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,MeshType::dim()> &xi);
};

/** Partial specialization for 2D meshes
 */
template<typename MeshType,typename BasisFunctionType>
class FunctionSpaceStructuredCheckNeighbouringElements<MeshType,BasisFunctionType,Mesh::isDim<2,MeshType>> :
  public FunctionSpaceStructuredFindPositionBase<MeshType,BasisFunctionType>
{
public:

  //! inherit constructor
  using FunctionSpaceStructuredFindPositionBase<MeshType,BasisFunctionType>::FunctionSpaceStructuredFindPositionBase;

protected:

  //! check if the point is in a neighbouring element to elementNo on ghostMeshNo (-1=main mesh, 0-5=ghost mesh on respective face, 0=face0Minus, 1=face0Plus, etc.), return true if the element was found amoung the neighbours
  //! set elementNo, ghostMeshNo and xi appropriately
  bool checkNeighbouringElements(const Vec3 &point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,MeshType::dim()> &xi);
};

/** Partial specialization for 3D meshes
 */
template<typename MeshType,typename BasisFunctionType>
class FunctionSpaceStructuredCheckNeighbouringElements<MeshType,BasisFunctionType,Mesh::isDim<3,MeshType>> :
  public FunctionSpaceStructuredFindPositionBase<MeshType,BasisFunctionType>
{
public:

  //! inherit constructor
  using FunctionSpaceStructuredFindPositionBase<MeshType,BasisFunctionType>::FunctionSpaceStructuredFindPositionBase;

protected:

  //! check if the point is in a neighbouring element to elementNo on ghostMeshNo (-1=main mesh, 0-5=ghost mesh on respective face, 0=face0Minus, 1=face0Plus, etc.), return true if the element was found amoung the neighbours
  //! set elementNo, ghostMeshNo and xi appropriately
  bool checkNeighbouringElements(const Vec3 &point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,MeshType::dim()> &xi);
};


}  // namespace

#include "function_space/09_function_space_structured_check_neighbouring_elements.tpp"
