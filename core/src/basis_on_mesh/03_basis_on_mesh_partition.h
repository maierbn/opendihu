#pragma once

#include <Python.h>  // has to be the first included header

#include "basis_on_mesh/02_basis_on_mesh_jacobian.h"
#include "partition/partition.h"

namespace BasisOnMesh
{
 
// forward declaration
template<typename MeshType, typename BasisFunctionType>
class BasisOnMesh {};
 
/** This adds functionality to create a partition / domain decomposition
 */
template<typename MeshType,typename BasisFunctionType>
class BasisOnMeshPartitionBase :
  public BasisOnMeshJacobian<MeshType,BasisFunctionType>
{
public:
  
  //! constructor
  BasisOnMeshPartition(std::shared_ptr<Partition::Manager> partitionManager);
  
  //! get the partition
  Partition::MeshPartition &meshPartition();
  
  //! get the partition as pointer of type meshPartitionBase, this is in the itnerface in mesh
  std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBase();

private:
  std::shared_ptr<Partition::Manager> partitionManager_;  ///< the partition manager object that can create partitions
  Partition::MeshPartition<BasisOnMesh<MeshType,BasisFunctionType>> meshPartition_;   ///< the partition information that is stored locally, i.e. the subdomain of the domain decomposition

};

/** specialization for structured meshes 
 */
template<typename MeshType,typename BasisFunctionType>
class BasisOnMeshPartition :
  public BasisOnMeshPartitionBase<MeshType,BasisFunctionType>
{
public:
  //! use inherited constructor 
  using BasisOnMeshPartitionBase<MeshType,BasisFunctionType>::BasisOnMeshPartitionBase;
  
  //! initiate the partitoning and then call the downwards initialize
  void initialize();

};

/** partial specialization for unstructured mesh
 */
template<int D,typename BasisFunctionType>
class BasisOnMeshGeometryData<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType> :
  public BasisOnMeshPartitionBase<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>
{
public:
  //! use inherited constructor
  using BasisOnMeshPartitionBase<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::BasisOnMeshPartitionBase;
  
  //! initiate the partitoning and then call the downwards initialize
  void initialize();

}  // namespace

#include "basis_on_mesh/03_basis_on_mesh_partition.tpp"
