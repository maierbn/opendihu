#include "function_space/10_function_space_field_variable.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"
#include "mesh/face_t.h"

namespace FunctionSpace
{

//! create a non-geometry field field variable with no values being set, with given component names
template<typename MeshType, typename BasisFunctionType>
std::shared_ptr<FieldVariable::FieldVariableBaseFunctionSpace<FunctionSpace<MeshType,BasisFunctionType>>> FunctionSpaceFieldVariable<MeshType,BasisFunctionType>::
createFieldVariable(std::string name, std::vector<std::string> componentNames)
{
  // create the field variable with template parameter nComponents by a factory class that perform the dynamic->static conversion
  std::shared_ptr<FieldVariable::FieldVariableBaseFunctionSpace<FunctionSpace<MeshType,BasisFunctionType>>> fieldVariable
    = FieldVariable::Factory<FunctionSpace<MeshType,BasisFunctionType>>::createFromFieldVariable(this->geometryField(), name, componentNames);

  // set entries to 0
  fieldVariable->zeroEntries();

  return fieldVariable;
}

//! create a non-geometry field field variable with no values being set, with given number of components, the component names will be the numbers
template<typename MeshType, typename BasisFunctionType>
std::shared_ptr<FieldVariable::FieldVariableBaseFunctionSpace<FunctionSpace<MeshType,BasisFunctionType>>> FunctionSpaceFieldVariable<MeshType,BasisFunctionType>::
createFieldVariable(std::string name, int nComponents)
{
  // create standard component names, the strings "0","1","2",...
  std::vector<std::string> componentNames(nComponents);
  for (int i=0; i<nComponents; i++)
  {
    componentNames[i] = std::to_string(i);
  }
  std::shared_ptr<FieldVariable::FieldVariableBaseFunctionSpace<FunctionSpace<MeshType,BasisFunctionType>>> fieldVariable
    = this->createFieldVariable(name, componentNames);

  // set entries to 0
  fieldVariable->zeroEntries();

  return fieldVariable;
}

//! create a non-geometry field field variable with no values being set, with given number of components, the component names will be the numbers
template<typename MeshType, typename BasisFunctionType>
template <int nComponents>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace<MeshType,BasisFunctionType>,nComponents>> FunctionSpaceFieldVariable<MeshType,BasisFunctionType>::
createFieldVariable(std::string name)
{
  // create standard component names, the strings "0","1","2",...
  std::vector<std::string> componentNames(nComponents);
  for (int i=0; i<nComponents; i++)
  {
    componentNames[i] = std::to_string(i);
  }
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace<MeshType,BasisFunctionType>,nComponents>> fieldVariable
    = std::make_shared<FieldVariable::FieldVariable<FunctionSpace<MeshType,BasisFunctionType>,nComponents>>(this->geometryField(), name, componentNames);

  // set entries to 0
  fieldVariable->zeroEntries();

  return fieldVariable;
}

//! create a non-geometry field field variable with no values being set, with given number of components and component names
template<typename MeshType, typename BasisFunctionType>
template <int nComponents>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace<MeshType,BasisFunctionType>,nComponents>> FunctionSpaceFieldVariable<MeshType,BasisFunctionType>::
createFieldVariable(std::string name, std::vector<std::string> componentNames)
{
  assert(this->meshPartition());
  assert(this->geometryField().functionSpace());
  assert(this->geometryField().functionSpace()->meshPartition());

  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace<MeshType,BasisFunctionType>,nComponents>> fieldVariable
    = std::make_shared<FieldVariable::FieldVariable<FunctionSpace<MeshType,BasisFunctionType>,nComponents>>(this->geometryField(), name, componentNames);

  // set entries to 0
  fieldVariable->zeroEntries();

  return fieldVariable;
}

//! create a non-geometry field field variable with no values being set, with given number of components, the component names will be the numbers
template<typename MeshType, typename BasisFunctionType>
int FunctionSpaceFieldVariable<MeshType,BasisFunctionType>::
getNumberScaleFactors(element_no_t elementGlobalNo)
{
  assert(this->geometryField_);
  return this->geometryField_->getNumberScaleFactors(elementGlobalNo);
}
  
template<typename MeshType, typename BasisFunctionType>
std::array<std::array<double,MeshType::dim()>,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> FunctionSpaceFieldVariable<MeshType,BasisFunctionType>::
getGradPhi(std::array<double,MeshType::dim()> xi) const
{
  // column-major storage, gradPhi[column][row] = gradPhi[dofIndex][i] = dphi_dofIndex/dxi_i, columnIdx = dofIndex, rowIdx = which direction
  std::array<std::array<double,MeshType::dim()>,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> gradPhi;
  for (int dofIndex = 0; dofIndex < FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement(); dofIndex++)
  {
    gradPhi[dofIndex] = this->gradPhi(dofIndex, xi);
  }
  return gradPhi;
}

template<typename MeshType, typename BasisFunctionType>
template <int nComponents, typename double_v_t>
std::array<double_v_t,nComponents> FunctionSpaceFieldVariable<MeshType,BasisFunctionType>::
interpolateValueInElement(std::array<std::array<double_v_t,nComponents>,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> &elementalDofValues,
                          std::array<double,MeshType::dim()> xi) const
{
  std::array<double_v_t,nComponents> result({0.0});
  for (int dofIndex = 0; dofIndex < FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement(); dofIndex++)
  {
    result += elementalDofValues[dofIndex]*this->phi(dofIndex,xi);
  }
  return result;
}

template<typename MeshType, typename BasisFunctionType>
template <typename double_v_t>
double_v_t FunctionSpaceFieldVariable<MeshType,BasisFunctionType>::
interpolateValueInElement(std::array<double_v_t,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> &elementalDofValues,
                          std::array<double,MeshType::dim()> xi) const
{
  double_v_t result{};
  for (int dofIndex = 0; dofIndex < FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement(); dofIndex++)
  {
    result += elementalDofValues[dofIndex]*this->phi(dofIndex,xi);
  }
  return result;
}

template<typename MeshType, typename BasisFunctionType>
std::array<double,MeshType::dim()> FunctionSpaceFieldVariable<MeshType,BasisFunctionType>::
interpolateGradientInElement(std::array<double,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> &elementalDofValues,
                             Tensor2<MeshType::dim()> inverseJacobianParameterSpace, std::array<double,MeshType::dim()> xi) const
{
  const int D = MeshType::dim();
  const int nDofsPerElement = FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement();

  std::array<double,D> gradPhiWorldSpace{0.0};

  // loop over dofs in element that contribute to the gradient at dofIndex
  for (int dofIndex = 0; dofIndex < nDofsPerElement; dofIndex++)
  {
    // get gradient at dof
    std::array<double,D> gradPhiParameterSpace = this->gradPhi(dofIndex, xi);

    //VLOG(2) << "  dofIndex2=" << dofIndex << ", xi=" << xi << ", gradPhiParameterSpace: " << gradPhiParameterSpace;


    std::array<double,D> gradPhiWorldSpaceDofIndex2{0.0};

    // transform grad from parameter space to world space
    for (int direction = 0; direction < D; direction++)
    {
      //VLOG(2) << "   component " << direction;
      for (int k = 0; k < D; k++)
      {
        // jacobianParameterSpace[columnIdx][rowIdx] = dX_rowIdx/dxi_columnIdx
        // inverseJacobianParameterSpace[columnIdx][rowIdx] = dxi_rowIdx/dX_columnIdx because of inverse function theorem

        const double dphiDofIndex2_dxik = gradPhiParameterSpace[k];   // dphi_dofIndex/dxi_k
        const double dxik_dXdirection = inverseJacobianParameterSpace[direction][k];  // dxi_k/dX_direction


        //VLOG(2) << "     += " << dphiDofIndex2_dxik << " * " << dxik_dXdirection;

        gradPhiWorldSpaceDofIndex2[direction] += dphiDofIndex2_dxik * dxik_dXdirection;
      }
    }

    //VLOG(2) << "  gradPhiWorldSpaceDofIndex2: " << gradPhiWorldSpaceDofIndex2
    //  << " multiply with solution value at dof " << dofIndex << ", " << elementalDofValues[dofIndex];

    //VLOG(2) << " sum contributions from the other ansatz functions at this dof: " << gradPhiWorldSpace;

    gradPhiWorldSpace += gradPhiWorldSpaceDofIndex2 * elementalDofValues[dofIndex];

    //VLOG(2) << "                                                             -> " << gradPhiWorldSpace;
  }  // dofIndex

  return gradPhiWorldSpace;
}

template<typename MeshType, typename BasisFunctionType>
template<typename double_v_t>
VecD<3,double_v_t> FunctionSpaceFieldVariable<MeshType,BasisFunctionType>::
getNormal(Mesh::face_t face, std::array<VecD<3,double_v_t>,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> geometryValues,
          std::array<double,MeshType::dim()> xi)
{
  // compute normal analog to nansons formula
  // Nansons formula: ds = J F^-T dS (ds, dS are normal vectors, here ds is in world space, dS is in index space)
  using Vec3v = VecD<3,double_v_t>;
  const int D = MeshType::dim();

  // compute the 3xD jacobian of the parameter space to world space mapping
  std::array<Vec3v,D> jacobian = this->computeJacobian(geometryValues, xi);
  std::array<Vec3v,3> jacobian3x3 = MathUtility::transformToDxD<3,D>(jacobian);

  // compute J F^-T, J = det F, F = jacobian
  std::array<Vec3v,3> cofactor = MathUtility::computeCofactorMatrix<double_v_t>(jacobian3x3);

  // transform the index space normal using Nanson's formula
  Vec3 normalIndexSpace = MathUtility::transformToD<3,D>(Mesh::getNormal<D>(face));
  Vec3v result = cofactor * normalIndexSpace;

  LOG(DEBUG) << "geometryValues: " << geometryValues;
  LOG(DEBUG) << "jacobian: " << jacobian << ", jacobian3x3: " << jacobian3x3;
  LOG(DEBUG) << "normal index space: " << Mesh::getNormal<D>(face) << ", " << normalIndexSpace;
  LOG(DEBUG) << "cofactor: " << cofactor << ", result: " << result;
  MathUtility::normalize<3>(result);
  return result;
}

template<typename MeshType, typename BasisFunctionType>
Vec3 FunctionSpaceFieldVariable<MeshType,BasisFunctionType>::
getNormal(Mesh::face_t face, element_no_t elementNoLocal, std::array<double,MeshType::dim()> xi)
{
  // compute normal analoguous to nansons formula
  // Nansons formula: ds = J F^-T dS (ds, dS are normal vectors, here ds is in world space, dS is in index space)

  // get geometry field values of element
  std::array<Vec3,FunctionSpaceBaseDim<MeshType::dim(),BasisFunctionType>::nDofsPerElement()> geometryValues;
  this->getElementGeometry(elementNoLocal, geometryValues);

  //LOG(DEBUG) << "elementNoLocal: " << elementNoLocal << ", geometryValues: " << geometryValues;
  return getNormal(face, geometryValues, xi);
}

template<typename MeshType, typename BasisFunctionType>
VecD<3,Vc::double_v> FunctionSpaceFieldVariable<MeshType,BasisFunctionType>::
getNormal(Mesh::face_t face, Vc::int_v elementNoLocal, std::array<double,MeshType::dim()> xi)
{
  // compute normal analoguous to nansons formula
  // Nansons formula: ds = J F^-T dS (ds, dS are normal vectors, here ds is in world space, dS is in index space)

  // get geometry field values of element
  std::array<VecD<3,Vc::double_v>,FunctionSpaceBaseDim<MeshType::dim(),BasisFunctionType>::nDofsPerElement()> geometryValues;
  this->getElementGeometry(elementNoLocal, geometryValues);

  //LOG(DEBUG) << "elementNoLocal: " << elementNoLocal << ", geometryValues: " << geometryValues;
  return getNormal(face, geometryValues, xi);
}

template<typename MeshType, typename BasisFunctionType>
template<typename double_v_t, typename element_no_v_t>
Tensor2<MeshType::dim(),double_v_t> FunctionSpaceFieldVariable<MeshType,BasisFunctionType>::
getInverseJacobian(std::array<VecD<3,double_v_t>,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()> &geometryValues, element_no_v_t elementNo, std::array<double,MeshType::dim()> xi)
{
  // define constants
  const int D = MeshType::dim();

  // compute the 3xD jacobian of the parameter space to world space mapping
  Tensor2<D,double_v_t> jacobianParameterSpace = MathUtility::transformToDxD<D,D>(this->computeJacobian(geometryValues, xi));
  double_v_t jacobianDeterminant;
  double_v_t approximateMeshWidth = MathUtility::computeApproximateMeshWidth<double_v_t,FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement()>(geometryValues);
  Tensor2<D,double_v_t> inverseJacobianParameterSpace = MathUtility::computeInverse(jacobianParameterSpace, approximateMeshWidth, jacobianDeterminant);

  return inverseJacobianParameterSpace;
}

template<typename MeshType, typename BasisFunctionType>
Tensor2<MeshType::dim()> FunctionSpaceFieldVariable<MeshType,BasisFunctionType>::
getInverseJacobian(element_no_t elementNo, std::array<double,MeshType::dim()> xi)
{
  const int nDofsPerElement = FunctionSpaceFunction<MeshType,BasisFunctionType>::nDofsPerElement();
  
  // get geometry field (which are the node positions for Lagrange basis and node positions and derivatives for Hermite)
  std::array<Vec3,nDofsPerElement> geometryValues;
  this->mesh_->getElementGeometry(elementNo, geometryValues);

  return getInverseJacobian(geometryValues, elementNo, xi);
}

} // namespace
