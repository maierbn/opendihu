#pragma once

#include <Python.h>  // has to be the first included header

#include "basis_on_mesh/07_basis_on_mesh_nodes.h"

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
                                 
  //! interpolate the value within an element at the given xi position using the basis functions
  double interpolateValueInElement(std::array<double,BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement()> &elementalDofValues,
                                   std::array<double,MeshType::dim()> xi) const;
                                    
  //! compute the normal in world space, normal to face at xi, use the given geometry values, that can by obtained by fieldVariable->getElementValues(elementNo, geometryValues) or mesh->getElementGeometry(elementNo, geometryValues)
  Vec3 getNormal(Mesh::face_t face, std::array<Vec3,BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement()> geometryValues, std::array<double,MeshType::dim()> xi);
  
  //! compute the normal in world space, normal to face at xi
  Vec3 getNormal(Mesh::face_t face, element_no_t elementNo, std::array<double,MeshType::dim()> xi);
  
  //! get the element no and the xi value of the point, return true if the point is inside the mesh or false otherwise
  bool findPosition(Vec3 point, element_no_t &elementNo, std::array<double,MeshType::dim()> &xi);
  
  //! check if the point lies inside the element, if yes, return true and set xi to the value of the point
  bool pointIsInElement(Vec3 point, element_no_t elementNo, std::array<double,MeshType::dim()> &xi);
};

}  // namespace

#include "basis_on_mesh/08_basis_on_mesh_field_variable.tpp"
