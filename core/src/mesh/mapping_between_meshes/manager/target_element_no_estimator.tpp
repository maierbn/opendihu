#include "mesh/mapping_between_meshes/manager/target_element_no_estimator.h"

#include <Python.h>  // has to be the first included header

namespace MappingBetweenMeshes
{

template<typename BasisFunctionType>
TargetElementNoEstimator<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunctionType>>::
TargetElementNoEstimator(std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunctionType>> sourceFunctionSpace,
                         std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunctionType>> targetFunctionSpace) :
  sourceFunctionSpace_(sourceFunctionSpace), targetFunctionSpace_(targetFunctionSpace), targetElementForYChange_(0), targetElementForZChange_(0)
{
}

//! try to improve the estimation for targetElementNo, in which sourceDofNoLocal should be
template<typename BasisFunctionType>
void TargetElementNoEstimator<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunctionType>>::
estimateElementNo(node_no_t sourceDofNoLocal, element_no_t &targetElementNo)
{
  // get local coordintes for sourceDofNoLocal
  node_no_t nodeNoLocal = sourceDofNoLocal / sourceFunctionSpace_->nDofsPerNode();
  std::array<int,3> coordinatesSource = sourceFunctionSpace_->meshPartition()->getCoordinatesLocal(nodeNoLocal);

  // sourceFunctionSpace_->meshPartition->nNodesLocalWithoutGhosts()-1

  // if source dof begins new x row
  if (coordinatesSource[0] == 0)
  {
    targetElementNo = targetElementForYChange_;
  }
  else if (coordinatesSource[0] == 1)
  {
    targetElementForYChange_ = targetElementNo;
  }

  // if source dof begins new x-y plane, i.e. z increases
  if (coordinatesSource[1] == 0 && coordinatesSource[0] == 0)
  {
    targetElementNo = targetElementForZChange_;
  }
  else if (coordinatesSource[1] == 0 && coordinatesSource[0] == 1)
  {
    targetElementForZChange_ = targetElementNo;
  }
}

template<typename BasisFunctionType>
TargetElementNoEstimator<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<3>, BasisFunctionType>>::
TargetElementNoEstimator(std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunctionType>> sourceFunctionSpace,
                         std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<3>, BasisFunctionType>> targetFunctionSpace) :
  sourceFunctionSpace_(sourceFunctionSpace), targetFunctionSpace_(targetFunctionSpace), targetElementForYChange_(0), targetElementForZChange_(0)
{
}


//! try to improve the estimation for targetElementNo, in which sourceDofNoLocal should be
template<typename BasisFunctionType>
void TargetElementNoEstimator<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunctionType>, FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<3>, BasisFunctionType>>::
estimateElementNo(node_no_t sourceDofNoLocal, element_no_t &targetElementNo)
{
  /*
  // get local coordintes for sourceDofNoLocal
  int subMeshNo = targetFunctionSpace_->subMeshNoWherePointWasFound();
  using SubMeshType = FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunctionType>;

  if (subMeshNo < 0 || subMeshNo >= targetFunctionSpace_->subFunctionSpaces().size())
    return;

  const std::shared_ptr<SubMeshType> &subFunctionSpaceTarget = targetFunctionSpace_->subFunctionSpaces()[subMeshNo];
  */

  // get local coordintes for sourceDofNoLocal
  node_no_t nodeNoLocal = sourceDofNoLocal / sourceFunctionSpace_->nDofsPerNode();
  std::array<int,3> coordinatesSource = sourceFunctionSpace_->meshPartition()->getCoordinatesLocal(nodeNoLocal);

  // sourceFunctionSpace_->meshPartition->nNodesLocalWithoutGhosts()-1

  // if source dof begins new x row
  if (coordinatesSource[0] == 0)
  {
    targetElementNo = targetElementForYChange_;
  }
  else if (coordinatesSource[0] == 1)
  {
    targetElementForYChange_ = targetElementNo;
  }

  // if source dof begins new x-y plane, i.e. z increases
  if (coordinatesSource[1] == 0 && coordinatesSource[0] == 0)
  {
    targetElementNo = targetElementForZChange_;
  }
  else if (coordinatesSource[1] == 0 && coordinatesSource[0] == 1)
  {
    targetElementForZChange_ = targetElementNo;
  }
}

}  // namespace
