#include "function_space/03_function_space_partition.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"
#include "utility/python_utility.h"

namespace FunctionSpace
{

// forward declaration
template<typename MeshType,typename BasisFunctionType>
class FunctionSpace;

template<int D,typename BasisFunctionType>
FunctionSpacePartition<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
FunctionSpacePartition(std::shared_ptr<Partition::Manager> partitionManager,
                       std::vector<std::shared_ptr<FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>> subFunctionSpaces) :
  FunctionSpacePartitionBase<Mesh::CompositeOfDimension<D>,BasisFunctionType>::FunctionSpacePartitionBase(partitionManager, PythonConfig()),
  subFunctionSpaces_(subFunctionSpaces)
{
  // the subFunctionSpaces have been created and initialized by the mesh manager
}

template<int D,typename BasisFunctionType>
void FunctionSpacePartition<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
initialize()
{
  VLOG(1) << "FunctionSpacePartition<Composite>::initialize()";

  // initialize number of local and global elements
  this->nElementsLocal_ = 0;
  this->nElementsGlobal_ = 0;

  // iterate over submeshes
  for(std::shared_ptr<SubFunctionSpaceType> &subFunctionSpace : subFunctionSpaces_)
  {
    subFunctionSpace->initialize();

    // count number of elements
    this->nElementsLocal_ += subFunctionSpace->nElementsLocal();
    this->nElementsGlobal_ += subFunctionSpace->nElementsGlobal();
  }

  // create partitioning
  assert(this->partitionManager_ != nullptr);
  std::vector<int> rankNos;
  if (this->specificSettings_.hasKey("rankNos"))
    this->specificSettings_.template getOptionVector<int>("rankNos", rankNos);

  this->meshPartition_ = this->partitionManager_->template createPartitioningComposite<BasisFunctionType,D>(subFunctionSpaces_, rankNos);

  assert(this->meshPartition_);

  // set initalized_ to true which indicates that initialize has been called
  this->initialized_ = true;
}

template<int D,typename BasisFunctionType>
const std::vector<std::shared_ptr<FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>> &
FunctionSpacePartition<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
subFunctionSpaces()
{
  return subFunctionSpaces_;
}

} // namespace
