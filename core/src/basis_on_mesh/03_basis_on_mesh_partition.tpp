#include "basis_on_mesh/03_basis_on_mesh_partition.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"

namespace BasisOnMesh
{
template<typename MeshType,typename BasisFunctionType>
BasisOnMeshPartition<MeshType,BasisFunctionType>:: 
BasisOnMeshPartition(std::shared_ptr<Partition::Manager> partitionManager, PyObject *specificSettings) : 
  partitionManager_(partitionManager), BasisOnMeshJacobian<MeshType,BasisFunctionType>(specificSettings)
{
}

template<typename MeshType,typename BasisFunctionType>
void BasisOnMeshPartition<MeshType,BasisFunctionType>::
initialize()
{
  // Creation of the partitioning is only possible after the number of elements is known.
  // Because this may need file I/O (e.g. reading from exfiles)
 
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