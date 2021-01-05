#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"
#include "function_space/08_function_space_elemental_dofs.h"
#include "function_space/00_function_space_base_dim.h"
#include "mesh/mesh.h"

namespace FunctionSpace
{

/** Adds functionality to get local nodes nos of the elements
 */
template<typename MeshType,typename BasisFunctionType>
class FunctionSpaceNodes :
  public FunctionSpaceElementalDofs<MeshType,BasisFunctionType>
{
public:

  //! inherit constructor
  using FunctionSpaceElementalDofs<MeshType,BasisFunctionType>::FunctionSpaceElementalDofs;

  //! return an array of all node nos. of the element
  std::array<dof_no_t,FunctionSpaceFunction<MeshType,BasisFunctionType>::nNodesPerElement()>
  getElementNodeNos(element_no_t elementNo) const;

  //! get the face that is defined by the dof nos in the element
  Mesh::face_t getFaceFromElementalDofNos(std::array<int,FunctionSpaceBaseDim<MeshType::dim()-1,BasisFunctionType>::nDofsPerElement()> elementalDofNos);
};

}  // namespace

#include "function_space/09_function_space_nodes.tpp"
