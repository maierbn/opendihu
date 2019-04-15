#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"
#include "function_space/07_function_space_faces.h"
#include "mesh/mesh.h"

namespace FunctionSpace
{

/** FunctionSpace derives from the mesh and adds functionality to get dof and node numbers, phi and gradient.
 */
template<typename MeshType,typename BasisFunctionType>
class FunctionSpaceNodes :
  public FunctionSpaceFaces<MeshType,BasisFunctionType>
{
public:

  //! inherit constructor
  using FunctionSpaceFaces<MeshType,BasisFunctionType>::FunctionSpaceFaces;

  //! return an array of all node nos. of the element
  std::array<dof_no_t,FunctionSpaceFunction<MeshType,BasisFunctionType>::nNodesPerElement()>
  getElementNodeNos(element_no_t elementNo) const;

};

/** Partial specialization for CompletePolynomials which do not need nodes and thus have no nodes functionality.
 */
template<typename MeshType,int D,int order>
class FunctionSpaceNodes<MeshType, BasisFunction::CompletePolynomialOfDimensionAndOrder<D,order>> :
  public FunctionSpaceFunction<MeshType, BasisFunction::CompletePolynomialOfDimensionAndOrder<MeshType::dim(),order>>
{
public:

  //! constructor, with second bool argument, inherits from FunctionSpaceFunction
  FunctionSpaceNodes(PythonConfig specificSettings, bool noGeometryField=false);

  //! get the number of dofs
  dof_no_t nDofsLocal() const;

};

}  // namespace

#include "function_space/08_function_space_nodes.tpp"
