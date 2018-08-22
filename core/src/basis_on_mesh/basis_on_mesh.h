#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"
#include "basis_on_mesh/10_basis_on_mesh_xi.h"
#include "mesh/mesh.h"

namespace BasisOnMesh
{

/** BasisOnMesh derives from the mesh and adds functionality to get dof and node numbers, phi and gradient.
 */
template<typename MeshType,typename BasisFunctionType>
class BasisOnMesh : public BasisOnMeshXi<MeshType,BasisFunctionType>
{
public:

  //! inherit constructor
  using BasisOnMeshXi<MeshType,BasisFunctionType>::BasisOnMeshXi;

  typedef MeshType Mesh;
  typedef BasisFunctionType BasisFunction;
  typedef BasisOnMesh<MeshType,BasisFunctionType> HighOrderBasisOnMesh;

  //! return an array of all dof nos. of the element (local dof nos)
  std::array<dof_no_t,BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement()>
  getElementDofNosLocal(element_no_t elementNo) const;

  //! set a vector of all dof nos. of the element
  void getElementDofNosLocal(element_no_t elementNo, std::vector<dof_no_t> &globalDofNos) const;
};

/** Partial specialization for CompletePolynomials which do not need nodes and thus have no nodes functionality.
 */
template<typename MeshType,int D,int order>
class BasisOnMesh<MeshType, BasisFunction::CompletePolynomialOfDimensionAndOrder<D,order>> :
  public BasisOnMeshXi<MeshType, BasisFunction::CompletePolynomialOfDimensionAndOrder<MeshType::dim(),order>>
{
public:

  //! inherit constructor
  using BasisOnMeshXi<MeshType, BasisFunction::CompletePolynomialOfDimensionAndOrder<MeshType::dim(),order>>::BasisOnMeshXi;

  typedef MeshType Mesh;
  typedef BasisFunction::CompletePolynomialOfDimensionAndOrder<D,order> BasisFunction;

  //! return an array of all dof nos. of the element
  std::array<dof_no_t,BasisOnMeshFunction<MeshType,BasisFunction>::nDofsPerElement()>
  getElementDofNosLocal(element_no_t elementNo) const;

  //! set a vector of all dof nos. of the element
  void getElementDofNosLocal(element_no_t elementNo, std::vector<dof_no_t> &globalDofNos) const;

  // the following methods are there for compatibility with the interface
  //! (unused method) initialize
  void initialize(){}

  //! (unused method) return number of nodes, 0 for this mesh
  node_no_t nNodesLocalWithGhosts() const {return 0;}

  //! (unused method) fill a vector with positions of the nodes, consecutive (x,y,z) values
  void getNodePositions(std::vector<double> &nodePositions) const {}

  //! (unused method) return the geometry field entry (node position for Lagrange elements) of a specific dof
  Vec3 getGeometry(node_no_t dofNo) const {return Vec3();}
};

template<typename BasisFunctionType>
class BasisOnMesh<Mesh::None, BasisFunctionType> : public Mesh::None
{
public:
  using None::None;
  void initialize(){}
};

}  // namespace

#include "basis_on_mesh/basis_on_mesh.tpp"
