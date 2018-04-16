#pragma once

#include <Python.h>  // has to be the first included header

#include "basis_on_mesh/06_basis_on_mesh_nodes.h"

namespace BasisOnMesh
{

//! forward declaration for field variables 
template<typename MeshType,typename BasisFunctionType>
class BasisOnMesh;

/** BasisOnMesh derives from the mesh and adds functionality to get dof and node numbers, phi and gradient.
 */
template<typename MeshType,typename BasisFunctionType>
class BasisOnMeshFieldVariable : 
  public BasisOnMeshNodes<MeshType,BasisFunctionType>
{
public:
   
  //! inherit constructor
  using BasisOnMeshNodes<MeshType,BasisFunctionType>::BasisOnMeshNodes;
  
  //! create a non-geometry field field variable with no values being set, with given component names
  std::shared_ptr<FieldVariable::FieldVariableBase<BasisOnMesh<MeshType,BasisFunctionType>>> createFieldVariable(std::string name, std::vector<std::string> componentNames);
  
  //! create a non-geometry field field variable with no values being set, with given number of components, the component names will be the numbers
  std::shared_ptr<FieldVariable::FieldVariableBase<BasisOnMesh<MeshType,BasisFunctionType>>> createFieldVariable(std::string name, int nComponents=1);
  
  //! create a non-geometry field field variable with no values being set, with given number of components, the component names will be the numbers
  template <int nComponents>
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMesh<MeshType,BasisFunctionType>,nComponents>> createFieldVariable(std::string name);
  
  //! return an array of the gradients of all nodal basis functions, evaluated at xi  
  std::array<std::array<double,MeshType::dim()>,BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement()> 
  getGradPhi(std::array<double,MeshType::dim()> xi) const;
  
  //! interpolate the nComponents values within an element at the given xi position using the basis functions
  template <int nComponents>
  std::array<double,nComponents> interpolateValueInElement(std::array<std::array<double,nComponents>,BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement()> &elementalDofValues,
                                                  std::array<double,MeshType::dim()> xi) const;
};

}  // namespace

#include "basis_on_mesh/07_basis_on_mesh_field_variable.tpp"
