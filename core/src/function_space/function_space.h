#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"
#include "mesh/type_traits.h"
#include "function_space/11_function_space_xi.h"
#include "mesh/mesh.h"
#include "mesh/surface_mesh.h"
#include "basis_function/lagrange.h"
#include "function_space/function_space_generic.h"
#include <vc_or_std_simd.h>  // this includes <Vc/Vc> or a Vc-emulating wrapper of <experimental/simd> if available

namespace FunctionSpace
{

/** FunctionSpace derives from the mesh and adds functionality to get dof and node numbers, phi and gradient.
 */
template<typename MeshType,typename BasisFunctionType>
class FunctionSpace : public FunctionSpacePointInElement<MeshType,BasisFunctionType>
{
public:

  //! inherit constructor
  using FunctionSpacePointInElement<MeshType,BasisFunctionType>::FunctionSpacePointInElement;

  typedef MeshType Mesh;
  typedef BasisFunctionType BasisFunction;
  typedef FunctionSpace<MeshType,BasisFunctionType> HighOrderFunctionSpace;
  typedef typename ::Mesh::SurfaceMesh<MeshType>::type SurfaceMesh;

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

  //! get a description of the function space, with mesh name and type
  std::string getDescription() const;
};

}  // namespace

#include "function_space/function_space.tpp"
