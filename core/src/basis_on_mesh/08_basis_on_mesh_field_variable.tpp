#include "basis_on_mesh/08_basis_on_mesh_field_variable.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"
#include "mesh/face_t.h"

namespace BasisOnMesh
{

//! create a non-geometry field field variable with no values being set, with given component names
template<typename MeshType, typename BasisFunctionType>
std::shared_ptr<FieldVariable::FieldVariableBase<BasisOnMesh<MeshType,BasisFunctionType>>> BasisOnMeshFieldVariable<MeshType,BasisFunctionType>::
createFieldVariable(std::string name, std::vector<std::string> componentNames)
{
  // create the field variable with template parameter nComponents by a factory class that perform the dynamic->static conversion
  std::shared_ptr<FieldVariable::FieldVariableBase<BasisOnMesh<MeshType,BasisFunctionType>>> fieldVariable
    = FieldVariable::Factory<BasisOnMesh<MeshType,BasisFunctionType>>::createFromFieldVariable(this->geometryField(), name, componentNames);
  
  return fieldVariable;
}

//! create a non-geometry field field variable with no values being set, with given number of components, the component names will be the numbers
template<typename MeshType, typename BasisFunctionType>
std::shared_ptr<FieldVariable::FieldVariableBase<BasisOnMesh<MeshType,BasisFunctionType>>> BasisOnMeshFieldVariable<MeshType,BasisFunctionType>::
createFieldVariable(std::string name, int nComponents)
{
  // create standard component names, the strings "0","1","2",...
  std::vector<std::string> componentNames(nComponents);
  for (int i=0; i<nComponents; i++)
  {
    componentNames[i] = std::to_string(i);
  }
  std::shared_ptr<FieldVariable::FieldVariableBase<BasisOnMesh<MeshType,BasisFunctionType>>> fieldVariable 
    = this->createFieldVariable(name, componentNames); 
  
  return fieldVariable;
}

//! create a non-geometry field field variable with no values being set, with given number of components, the component names will be the numbers
template<typename MeshType, typename BasisFunctionType>
template <int nComponents>
std::shared_ptr<FieldVariable::FieldVariable<BasisOnMesh<MeshType,BasisFunctionType>,nComponents>> BasisOnMeshFieldVariable<MeshType,BasisFunctionType>::
createFieldVariable(std::string name)
{
  // create standard component names, the strings "0","1","2",...
  std::vector<std::string> componentNames(nComponents);
  for (int i=0; i<nComponents; i++)
  {
    componentNames[i] = std::to_string(i);
  }
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMesh<MeshType,BasisFunctionType>,nComponents>> fieldVariable
    = std::make_shared<FieldVariable::FieldVariable<BasisOnMesh<MeshType,BasisFunctionType>,nComponents>>();
  
  fieldVariable->initializeFromFieldVariable(this->geometryField(), name, componentNames);
  return fieldVariable;
}
 
template<typename MeshType, typename BasisFunctionType>
std::array<std::array<double,MeshType::dim()>,BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement()> BasisOnMeshFieldVariable<MeshType,BasisFunctionType>::
getGradPhi(std::array<double,MeshType::dim()> xi) const
{
  // column-major storage, gradPhi[column][row] = gradPhi[dofIndex][i] = dphi_dofIndex/dxi_i, columnIdx = dofIndex, rowIdx = which direction
  std::array<std::array<double,MeshType::dim()>,BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement()> gradPhi;
  for (int dofIndex = 0; dofIndex < BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement(); dofIndex++)
  {
    gradPhi[dofIndex] = this->gradPhi(dofIndex, xi);
  }
  return gradPhi;
}

template<typename MeshType, typename BasisFunctionType>
template <int nComponents>
std::array<double,nComponents> BasisOnMeshFieldVariable<MeshType,BasisFunctionType>::
interpolateValueInElement(std::array<std::array<double,nComponents>,BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement()> &elementalDofValues,
                 std::array<double,MeshType::dim()> xi) const
{
  std::array<double,nComponents> result;
  for (int dofIndex = 0; dofIndex < BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement(); dofIndex++)
  {
    result += elementalDofValues[dofIndex]*this->phi(dofIndex,xi);
  }
  return result;
}

template<typename MeshType, typename BasisFunctionType>
Vec3 BasisOnMeshFieldVariable<MeshType,BasisFunctionType>::
getNormal(Mesh::face_t face, std::array<Vec3,BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement()> geometryValues, std::array<double,MeshType::dim()> xi)
{
  // compute normal analog to nansons formula
  // Nansons formula: ds = J F^-T dS (ds, dS are normal vectors, here ds is in world space, dS is in index space)
 
  const int D = MeshType::dim(); 
 
  // compute the 3xD jacobian of the parameter space to world space mapping
  std::array<Vec3,D> jacobian = this->computeJacobian(geometryValues, xi);
  
  std::array<Vec3,3> jacobian3x3 = MathUtility::transformToDxD<3,D>(jacobian);
      
  // compute J F^-T, J = det F, F = jacobian
  std::array<Vec3,3> cofactor = MathUtility::computeCofactorMatrix<3>(jacobian3x3);
  
  // transform the index space normal using Nanson's formula
  Vec3 normalIndexSpace = MathUtility::transformToD<3,D>(Mesh::getNormal<D>(face));
  Vec3 result = cofactor * normalIndexSpace;
  MathUtility::normalize<3>(result);
  return result;
}
  

template<typename MeshType, typename BasisFunctionType>
Vec3 BasisOnMeshFieldVariable<MeshType,BasisFunctionType>::
getNormal(Mesh::face_t face, element_no_t elementNo, std::array<double,MeshType::dim()> xi)
{
  // compute normal analog to nansons formula
  // Nansons formula: ds = J F^-T dS (ds, dS are normal vectors, here ds is in world space, dS is in index space)
 
  // get geometry field values of element
  std::array<Vec3,BasisOnMeshBaseDim<MeshType::dim(),BasisFunctionType>::nDofsPerElement()> geometryValues;
  this->getElementGeometry(elementNo, geometryValues);
  
  return getNormal(face, geometryValues, xi);
}

};  // namespace