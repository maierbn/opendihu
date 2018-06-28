#include "basis_on_mesh/03_basis_on_mesh_partition.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"

namespace BasisOnMesh
{
template<typename MeshType,typename BasisFunctionType>
BasisOnMeshPartition<MeshType,BasisFunctionType>:: 
BasisOnMeshPartition(std::shared_ptr<Partition::Manager> partitionManager) : 
  partitionManager_(partitionManager)
{
}
 
template<typename MeshType,typename BasisFunctionType>
void BasisOnMeshPartition<MeshType,BasisFunctionType>::
initialize()
{
  // create partitioning
  assert(partitionManager_ != nullptr);
  partition_ = partitionManager_->createPartitioning(this->nElements());
}

template<typename MeshType,typename BasisFunctionTypes>
Partition::MeshPartition &BasisOnMeshPartition<MeshType,BasisFunctionType>::
partition()
{
  return partition_;
}


};  // namespace