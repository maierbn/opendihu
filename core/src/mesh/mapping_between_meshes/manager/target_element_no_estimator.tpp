#include "mesh/mapping_between_meshes/manager/target_element_no_estimator.h"

#include <Python.h>  // has to be the first included header

namespace MappingBetweenMeshes
{

template<typename BasisFunctionType1, typename BasisFunctionType2>
TargetElementNoEstimator<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunctionType1>, FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunctionType2>>::
TargetElementNoEstimator(std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunctionType1>> sourceFunctionSpace,
                         std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunctionType2>> targetFunctionSpace) :
  sourceFunctionSpace_(sourceFunctionSpace), targetFunctionSpace_(targetFunctionSpace), targetElementForYChange_(0), targetElementForZChange_(0)
{
}

//! try to improve the estimation for targetElementNo, in which sourceDofNoLocal should be
template<typename BasisFunctionType1, typename BasisFunctionType2>
void TargetElementNoEstimator<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunctionType1>, FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunctionType2>>::
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

template<typename BasisFunctionType1, typename BasisFunctionType2>
TargetElementNoEstimator<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunctionType1>, FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<3>, BasisFunctionType2>>::
TargetElementNoEstimator(std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunctionType1>> sourceFunctionSpace,
                         std::shared_ptr<FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<3>, BasisFunctionType2>> targetFunctionSpace) :
  sourceFunctionSpace_(sourceFunctionSpace), targetFunctionSpace_(targetFunctionSpace), targetElementForYChange_(0), targetElementForZChange_(0)
{
}


//! try to improve the estimation for targetElementNo, in which sourceDofNoLocal should be
template<typename BasisFunctionType1, typename BasisFunctionType2>
void TargetElementNoEstimator<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, BasisFunctionType1>, FunctionSpace::FunctionSpace<Mesh::CompositeOfDimension<3>, BasisFunctionType2>>::
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

  element_no_t lastElementNo = targetElementNo;

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

  // debugging output
  VLOG(1) << lastElementNo << "; " << coordinatesSource << "/("
    << sourceFunctionSpace_->meshPartition()->nNodesLocalWithoutGhosts(0) << ","
    << sourceFunctionSpace_->meshPartition()->nNodesLocalWithoutGhosts(1) << ","
    << sourceFunctionSpace_->meshPartition()->nNodesLocalWithoutGhosts(2) << ") "
    << "[y:" << targetElementForYChange_ << ",z:" << targetElementForZChange_ << "] -> " << targetElementNo;
}

}  // namespace
