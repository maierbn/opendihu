#include "basis_on_mesh/03_basis_on_mesh_partition.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"

namespace BasisOnMesh
{
template<typename MeshType,typename BasisFunctionType>
BasisOnMeshPartition<MeshType,BasisFunctionType>:: 
BasisOnMeshPartitionBase(std::shared_ptr<Partition::Manager> partitionManager, PyObject *specificSettings) : 
  partitionManager_(partitionManager), BasisOnMeshJacobian<MeshType,BasisFunctionType>(specificSettings)
{
}

template<typename MeshType,typename BasisFunctionTypes>
Partition::MeshPartition<BasisOnMeshType<MeshType,BasisFunctionType>> &BasisOnMeshPartitionBase<MeshType,BasisFunctionType>::
meshPartition()
{
  return meshPartition_;
}

template<typename MeshType,typename BasisFunctionTypes>
std::shared_ptr<Partition::MeshPartitionBase> BasisOnMeshPartitionBase<MeshType,BasisFunctionType>::
meshPartitionBase()
{
  return std::static_pointer_cast<Partition::MeshPartitionBase>(meshPartition_);
}


};  // namespace