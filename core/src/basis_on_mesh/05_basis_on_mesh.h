#pragma once

#include <Python.h>  // has to be the first included header

#include <array>
#include "control/types.h"
#include "basis_on_mesh/04_basis_on_mesh_nodes.h"
#include "mesh/mesh.h"

namespace BasisOnMesh
{

/** BasisOnMesh derives from the mesh and adds functionality to get dof and node numbers, phi and gradient.
 */
template<typename MeshType,typename BasisFunctionType>
class BasisOnMesh : public BasisOnMeshNodes<MeshType,BasisFunctionType>
{
public:
   
  //! inherit constructor
  using BasisOnMeshNodes<MeshType,BasisFunctionType>::BasisOnMeshNodes;
  
  typedef MeshType Mesh;
  typedef BasisFunctionType BasisFunction;
  
  //! return an array of all dof nos. of the element  
  std::array<int,BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement()> 
  getElementDofNos(element_no_t elementNo) const;

  //! return an array of all node nos. of the element  
  std::array<int,BasisOnMeshFunction<MeshType,BasisFunctionType>::nNodesPerElement()> 
  getElementNodeNos(element_no_t elementNo) const;
  
  //! return an array of the gradients of all nodal basis functions, evaluated at xi  
  std::array<std::array<double,MeshType::dim()>,BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement()> 
  getGradPhi(std::array<double,MeshType::dim()> xi) const;
  
  //! interpolate the nComponents values within an element at the given xi position using the basis functions
  template <int nComponents>
  double interpolateValue(std::array<std::array<double,nComponents>,BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement()> &nodalValues,
                          std::array<double,MeshType::dim()> xi) const; 
  
  //! create a non-geometry field field variable with no values being set, with given component names
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMesh<MeshType,BasisFunctionType>>> createFieldVariable(std::string name, std::vector<std::string> componentNames);
  
  //! create a non-geometry field field variable with no values being set, with given number of components, the component names will be the numbers
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMesh<MeshType,BasisFunctionType>>> createFieldVariable(std::string name, int nComponents=1);
  
};

template<typename BasisFunctionType>
class BasisOnMesh<Mesh::None, BasisFunctionType> : public Mesh::None
{
public:
  using None::None;
  void initialize(){}
};

}  // namespace

#include "basis_on_mesh/05_basis_on_mesh.tpp"
