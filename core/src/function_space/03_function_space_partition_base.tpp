#include "function_space/03_function_space_partition.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"

namespace FunctionSpace
{

template<typename MeshType,typename BasisFunctionType>
FunctionSpacePartitionBase<MeshType,BasisFunctionType>:: 
FunctionSpacePartitionBase(std::shared_ptr<Partition::Manager> partitionManager, PyObject *specificSettings) : 
  FunctionSpaceJacobian<MeshType,BasisFunctionType>(specificSettings), partitionManager_(partitionManager)
{
}

template<typename MeshType,typename BasisFunctionType>
std::shared_ptr<Partition::MeshPartition<FunctionSpace<MeshType,BasisFunctionType>,MeshType>> 
FunctionSpacePartitionBase<MeshType,BasisFunctionType>::
meshPartition() const
{
  return meshPartition_;
}

template<typename MeshType,typename BasisFunctionType>
std::shared_ptr<Partition::MeshPartitionBase> FunctionSpacePartitionBase<MeshType,BasisFunctionType>::
meshPartitionBase()
{
  return std::static_pointer_cast<Partition::MeshPartitionBase>(meshPartition_);
}


};  // namespace
