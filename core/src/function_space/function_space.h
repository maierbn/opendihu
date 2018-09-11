#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"
#include "function_space/10_function_space_xi.h"
#include "mesh/mesh.h"
#include "basis_function/lagrange.h"

namespace FunctionSpace
{

/** FunctionSpace derives from the mesh and adds functionality to get dof and node numbers, phi and gradient.
 */
template<typename MeshType,typename BasisFunctionType>
class FunctionSpace : public FunctionSpaceXi<MeshType,BasisFunctionType>
{
public:

  //! inherit constructor
  using FunctionSpaceXi<MeshType,BasisFunctionType>::FunctionSpaceXi;

  typedef MeshType Mesh;
  typedef BasisFunctionType BasisFunction;
  typedef FunctionSpace<MeshType,BasisFunctionType> HighOrderFunctionSpace;

  //! return an array of all dof nos. of the element (local dof nos)
  std::array<dof_no_t,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()>
  getElementDofNosLocal(element_no_t elementNo) const;

  //! set a vector of all dof nos. of the element
  void getElementDofNosLocal(element_no_t elementNo, std::vector<dof_no_t> &globalDofNos) const;
};

/** Partial specialization for CompletePolynomials which do not need nodes and thus have no nodes functionality.
 */
template<typename MeshType,int D,int order>
class FunctionSpace<MeshType, BasisFunction::CompletePolynomialOfDimensionAndOrder<D,order>> :
  public FunctionSpaceXi<MeshType, BasisFunction::CompletePolynomialOfDimensionAndOrder<MeshType::dim(),order>>
{
public:

  //! inherit constructor
  using FunctionSpaceXi<MeshType, BasisFunction::CompletePolynomialOfDimensionAndOrder<MeshType::dim(),order>>::FunctionSpaceXi;

  typedef MeshType Mesh;
  typedef BasisFunction::CompletePolynomialOfDimensionAndOrder<D,order> BasisFunction;

  //! return an array of all dof nos. of the element
  std::array<dof_no_t,FunctionSpaceFunction<MeshType,BasisFunction>::nDofsPerElement()>
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

// define generic function space without logical real world mesh presententation, that can be used for generic field variables.
// For example for MOR the reduced vectors do not live on any mesh, but they need a function space to be defined and such that output writers work.
typedef FunctionSpace<Mesh::StructuredRegularFixedOfDimension<1>,BasisFunction::LagrangeOfOrder<1>> Generic;

}  // namespace

#include "function_space/function_space.tpp"
