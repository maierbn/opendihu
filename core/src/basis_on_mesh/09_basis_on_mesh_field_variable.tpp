#include "basis_on_mesh/09_basis_on_mesh_field_variable.h"

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
    = std::make_shared<FieldVariable::FieldVariable<BasisOnMesh<MeshType,BasisFunctionType>,nComponents>>(this->geometryField(), name, componentNames);

  return fieldVariable;
}

//! create a non-geometry field field variable with no values being set, with given number of components and component names
template<typename MeshType, typename BasisFunctionType>
template <int nComponents>
std::shared_ptr<FieldVariable::FieldVariable<BasisOnMesh<MeshType,BasisFunctionType>,nComponents>> BasisOnMeshFieldVariable<MeshType,BasisFunctionType>::
createFieldVariable(std::string name, std::vector<std::string> componentNames)
{
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMesh<MeshType,BasisFunctionType>,nComponents>> fieldVariable
    = std::make_shared<FieldVariable::FieldVariable<BasisOnMesh<MeshType,BasisFunctionType>,nComponents>>();

  fieldVariable->initializeFromFieldVariable(this->geometryField(), name, componentNames);
  return fieldVariable;
}

//! create a non-geometry field field variable with no values being set, with given number of components, the component names will be the numbers
template<typename MeshType, typename BasisFunctionType>
int BasisOnMeshFieldVariable<MeshType,BasisFunctionType>::
getNumberScaleFactors(element_no_t elementGlobalNo)
{
  assert(this->geometryField_);
  return this->geometryField_->getNumberScaleFactors(elementGlobalNo);
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
  std::array<double,nComponents> result({0.0});
  for (int dofIndex = 0; dofIndex < BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement(); dofIndex++)
  {
    result += elementalDofValues[dofIndex]*this->phi(dofIndex,xi);
  }
  return result;
}

template<typename MeshType, typename BasisFunctionType>
double BasisOnMeshFieldVariable<MeshType,BasisFunctionType>::
interpolateValueInElement(std::array<double,BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement()> &elementalDofValues,
                 std::array<double,MeshType::dim()> xi) const
{
  double result = 0;
  for (int dofIndex = 0; dofIndex < BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement(); dofIndex++)
  {
    result += elementalDofValues[dofIndex]*this->phi(dofIndex,xi);
  }
  return result;
}

template<typename MeshType, typename BasisFunctionType>
std::array<double,MeshType::dim()> BasisOnMeshFieldVariable<MeshType,BasisFunctionType>::
interpolateGradientInElement(std::array<double,BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement()> &elementalDofValues,
                             Tensor2<MeshType::dim()> inverseJacobianParameterSpace, std::array<double,MeshType::dim()> xi) const
{
  const int D = MeshType::dim();
  const int nDofsPerElement = BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement();

  std::array<double,D> gradPhiWorldSpace{0.0};

  // loop over dofs in element that contribute to the gradient at dofIndex
  for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
  {
    // get gradient at dof
    std::array<double,D> gradPhiParameterSpace = this->gradPhi(dofIndex, xi);

    VLOG(2) << "  dofIndex2=" << dofIndex << ", xi=" << xi << ", gradPhiParameterSpace: " << gradPhiParameterSpace;


    std::array<double,D> gradPhiWorldSpaceDofIndex2{0.0};

    // transform grad from parameter space to world space
    for (int direction = 0; direction < D; direction++)
    {
      VLOG(2) << "   component " << direction;
      for (int k = 0; k < D; k++)
      {
        // jacobianParameterSpace[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx
        // inverseJacobianParameterSpace[columnIdx][rowIdx] = dxi_rowIdx/dX_columnIdx because of inverse function theorem

        const double dphiDofIndex2_dxik = gradPhiParameterSpace[k];   // dphi_dofIndex/dxi_k
        const double dxik_dXdirection = inverseJacobianParameterSpace[direction][k];  // dxi_k/dX_direction


        VLOG(2) << "     += " << dphiDofIndex2_dxik << " * " << dxik_dXdirection;

        gradPhiWorldSpaceDofIndex2[direction] += dphiDofIndex2_dxik * dxik_dXdirection;
      }
    }

    VLOG(2) << "  gradPhiWorldSpaceDofIndex2: " << gradPhiWorldSpaceDofIndex2 
      << " multiply with solution value at dof " << dofIndex << ", " << elementalDofValues[dofIndex];

    VLOG(2) << " sum contributions from the other ansatz functions at this dof: " << gradPhiWorldSpace;

    gradPhiWorldSpace += gradPhiWorldSpaceDofIndex2 * elementalDofValues[dofIndex];

    VLOG(2) << "                                                             -> " << gradPhiWorldSpace;
  }  // dofIndex

  return gradPhiWorldSpace;
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

template<typename MeshType, typename BasisFunctionType>
bool BasisOnMeshFieldVariable<MeshType,BasisFunctionType>::
findPosition(Vec3 point, element_no_t &elementNo, std::array<double,MeshType::dim()> &xi)
{
  const element_no_t nElements = this->nElementsLocal();
 
  // set starting no to 0 if it was not given and is thus arbitrarily initialized
  if (elementNo < 0 || elementNo >= nElements)
    elementNo = 0;
  
  // look in every element, starting at elementNo-2
  element_no_t elementNoStart = (elementNo - 2 + nElements) % nElements;
  element_no_t elementNoEnd = (elementNo - 3 + nElements) % nElements;
  
  VLOG(3) << "elementNoStart: " << elementNoStart << ", elementNoEnd: " << elementNoEnd << ", nElements: " << nElements;
  
  for (element_no_t currentElementNo = elementNoStart; currentElementNo != elementNoEnd; currentElementNo++)
  {
    if (currentElementNo == nElements)
    {
      currentElementNo = 0;
      if (elementNoEnd == currentElementNo)
        break;
    }
   
    
    if (this->pointIsInElement(point, currentElementNo, xi))
    {
      elementNo = currentElementNo;
      return true;
    }
  }
  return false;
}

template<typename MeshType, typename BasisFunctionType>
Tensor2<MeshType::dim()> BasisOnMeshFieldVariable<MeshType,BasisFunctionType>::
getInverseJacobian(std::array<Vec3,BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement()> &geometryValues, element_no_t elementNo, std::array<double,MeshType::dim()> xi)
{
  // define constants
  const int D = MeshType::dim();

  // compute the 3xD jacobian of the parameter space to world space mapping
  Tensor2<D> jacobianParameterSpace = MathUtility::transformToDxD<D,D>(this->computeJacobian(geometryValues, xi));
  double jacobianDeterminant;
  Tensor2<D> inverseJacobianParameterSpace = MathUtility::computeInverse<D>(jacobianParameterSpace, jacobianDeterminant);

  return inverseJacobianParameterSpace;
}

template<typename MeshType, typename BasisFunctionType>
Tensor2<MeshType::dim()> BasisOnMeshFieldVariable<MeshType,BasisFunctionType>::
getInverseJacobian(element_no_t elementNo, std::array<double,MeshType::dim()> xi)
{
  const int nDofsPerElement = BasisOnMeshFunction<MeshType,BasisFunctionType>::nDofsPerElement();
  
  // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
  std::array<Vec3,nDofsPerElement> geometryValues;
  this->mesh_->getElementGeometry(elementNo, geometryValues);

  return getInverseJacobian(geometryValues, elementNo, xi);
}

};  // namespace