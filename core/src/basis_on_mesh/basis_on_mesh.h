#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"
#include "basis_on_mesh/08_basis_on_mesh_field_variable.h"
#include "mesh/mesh.h"

namespace BasisOnMesh
{

/** BasisOnMesh derives from the mesh and adds functionality to get dof and node numbers, phi and gradient.
 */
template<typename MeshType,typename BasisFunctionType>
class BasisOnMesh : public BasisOnMeshFieldVariable<MeshType,BasisFunctionType>
{
public:
   
  //! inherit constructor
  using BasisOnMeshFieldVariable<MeshType,BasisFunctionType>::BasisOnMeshFieldVariable;
  
  typedef MeshType Mesh;
  typedef BasisFunctionType BasisFunction;
  
  //! return an array of all dof nos. of the element  
  std::array<dof_no_t,BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement()> 
  getElementDofNos(element_no_t elementNo) const;
  
  //! set a vector of all dof nos. of the element  
  void getElementDofNos(element_no_t elementNo, std::vector<dof_no_t> &globalDofNos) const;
};

/** Partial specialization for CompletePolynomials which do not need nodes and thus have no nodes functionality.
 */
template<typename MeshType,int D,int order>
class BasisOnMesh<MeshType, BasisFunction::CompletePolynomialOfDimensionAndOrder<D,order>> :
  public BasisOnMeshFieldVariable<MeshType, BasisFunction::CompletePolynomialOfDimensionAndOrder<MeshType::dim(),order>>
{
public:
   
  //! inherit constructor
  using BasisOnMeshFieldVariable<MeshType, BasisFunction::CompletePolynomialOfDimensionAndOrder<MeshType::dim(),order>>::BasisOnMeshFieldVariable;
  
  typedef MeshType Mesh;
  typedef BasisFunction::CompletePolynomialOfDimensionAndOrder<D,order> BasisFunction;
  
  //! return an array of all dof nos. of the element  
  std::array<dof_no_t,BasisOnMeshFunction<MeshType,BasisFunction>::nDofsPerElement()> 
  getElementDofNos(element_no_t elementNo) const;
  
  //! set a vector of all dof nos. of the element  
  void getElementDofNos(element_no_t elementNo, std::vector<dof_no_t> &globalDofNos) const;
  
  // the following methods are there for compatibility with the interface
  //! (unused method) initialize
  void initialize(){}
  
  //! (unused method) return number of nodes, 0 for this mesh
  node_no_t nNodes() const {return 0;}
  
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
