#include "basis_on_mesh/08_basis_on_mesh_partition.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"

namespace BasisOnMesh
{

template<typename MeshType,typename BasisFunctionType>
void BasisOnMeshPartition<MeshType,BasisFunctionType>::
initialize(std::shared_ptr<Partition::Manager> partitionManager)
{
  setupPartitioning();
  BasisOnMeshGeometry<MeshType,BasisFunctionType>::initialize();
}


template<typename MeshType,typename BasisFunctionType>
void BasisOnMeshPartition<MeshType,BasisFunctionType>::
setupPartitioning(std::shared_ptr<Partition::Manager> partitionManager)
{
  // create partitioning
  partition_ = partitionManager->createPartitioning(this->nElements());
}

template<typename MeshType,typename BasisFunctionTypes>
Partition::MeshPartition &BasisOnMeshPartition<MeshType,BasisFunctionType>::
partition()
{
  return partition_;
}
};  // namespace