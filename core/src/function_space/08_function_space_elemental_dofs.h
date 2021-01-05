#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"
#include "function_space/07_function_space_faces.h"
#include "mesh/mesh.h"

namespace FunctionSpace
{

/** Adds functionality to get all local dof nos of an element
 */
template<typename MeshType,typename BasisFunctionType>
class FunctionSpaceElementalDofs:
  public FunctionSpaceFaces<MeshType,BasisFunctionType>
{
public:

  //! inherit constructor
  using FunctionSpaceFaces<MeshType,BasisFunctionType>::FunctionSpaceFaces;

  //! return an array of all dof nos. of the element, including ghost dofs (local dof nos)
  std::array<dof_no_t,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()>
  getElementDofNosLocal(element_no_t elementNo) const;

  //! vectorized version of getElementDofNosLocal
  std::array<Vc::int_v,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()>
  getElementDofNosLocal(Vc::int_v elementNo) const;

  //! fill a vector of all local dof nos. of the element, without ghost dofs
  void getElementDofNosLocalWithoutGhosts(element_no_t elementNo, std::vector<dof_no_t> &dofNosLocal) const;

  //! fill a vector of all local dof nos. of the element, including ghost dofs
  void getElementDofNosLocal(element_no_t elementNo, std::vector<dof_no_t> &localDofNos) const;
};

}  // namespace

#include "function_space/08_function_space_elemental_dofs.tpp"
