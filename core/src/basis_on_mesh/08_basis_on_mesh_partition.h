#pragma once

#include <Python.h>  // has to be the first included header

#include "basis_on_mesh/07_basis_on_mesh_nodes.h"
#include "partition/partition.h"

namespace BasisOnMesh
{

/** This adds functionality to create a partition / domain decomposition
 */
template<typename MeshType,typename BasisFunctionType>
class BasisOnMeshPartition :
  public BasisOnMeshNodes<MeshType,BasisFunctionType>
{
public:

  //! inherit constructor
  using BasisOnMeshNodes<MeshType,BasisFunctionType>::BasisOnMeshNodes;

  //! get the partition
  Partition::MeshPartition &partition();
  
  //! initiate the partitoning and then call the downwards initialize
  void initialize();

private:
  //! create a partition
  void setupPartitioning(std::shared_ptr<Partition::Manager> partitionManager);
  
  Partition::MeshPartition partition_;   ///< the partition information that is stored locally, i.e. the subdomain of the domain decomposition


};

}  // namespace

#include "basis_on_mesh/08_basis_on_mesh_partition.tpp"
