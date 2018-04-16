#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"
#include "basis_on_mesh/05_basis_on_mesh_dofs_nodes.h"
#include "mesh/mesh.h"

namespace BasisOnMesh
{

/** BasisOnMesh derives from the mesh and adds functionality to get dof and node numbers, phi and gradient.
 */
template<typename MeshType,typename BasisFunctionType>
class BasisOnMeshNodes : 
  public BasisOnMeshDofsNodes<MeshType,BasisFunctionType>
{
public:
   
  //! inherit constructor
  using BasisOnMeshDofsNodes<MeshType,BasisFunctionType>::BasisOnMeshDofsNodes;
  
  //! return an array of all node nos. of the element  
  std::array<dof_no_t,BasisOnMeshFunction<MeshType,BasisFunctionType>::nNodesPerElement()> 
  getElementNodeNos(element_no_t elementNo) const;
  
};

/** Partial specialization for CompletePolynomials which do not need nodes and thus have no nodes functionality.
 */
template<typename MeshType,int D,int order>
class BasisOnMeshNodes<MeshType, BasisFunction::CompletePolynomialOfDimensionAndOrder<D,order>> :
  public BasisOnMeshFunction<MeshType, BasisFunction::CompletePolynomialOfDimensionAndOrder<MeshType::dim(),order>>
{
public:
   
  //! constructor, with second bool argument, inherits from BasisOnMeshFunction
  BasisOnMeshNodes(PyObject *specificSettings, bool noGeometryField=false);
  
  //! get the number of dofs
  dof_no_t nDofs() const;
  
};

}  // namespace

#include "basis_on_mesh/06_basis_on_mesh_nodes.tpp"
