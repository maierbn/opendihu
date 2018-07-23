#pragma once

#include <Python.h>  // has to be the first included header

#include "basis_on_mesh/02_basis_on_mesh_jacobian.h"
#include "partition/partition.h"

namespace BasisOnMesh
{
 
/** This adds functionality to create a partition / domain decomposition
 */
template<typename MeshType,typename BasisFunctionType>
class BasisOnMeshPartition :
  public BasisOnMeshJacobian<MeshType,BasisFunctionType>
{
public:
  
  //! constructor
  BasisOnMeshPartition(std::shared_ptr<Partition::Manager> partitionManager);
  
  //! initiate the partitoning and then call the downwards initialize
  void initialize();

  //! get the partition
  Partition::MeshPartition &meshPartition();
  
  //! get the partition as pointer of type meshPartitionBase, this is in the itnerface in mesh
  std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBase();

private:
  std::shared_ptr<Partition::Manager> partitionManager_;  ///< the partition manager object that can create partitions
  Partition::MeshPartition meshPartition_;   ///< the partition information that is stored locally, i.e. the subdomain of the domain decomposition

};

}  // namespace

#include "basis_on_mesh/03_basis_on_mesh_partition.tpp"
