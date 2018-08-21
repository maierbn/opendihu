#pragma once

#include <Python.h>  // has to be the first included header

#include "basis_on_mesh/02_basis_on_mesh_jacobian.h"
#include "partition/01_mesh_partition.h"
#include "partition/00_mesh_partition_base.h"

// forward declaration
namespace Partition 
{
class Manager;
//template<typename BasisOnMeshType, typename DummyForTraits>
//class MeshPartition;
}

namespace BasisOnMesh
{
 
// forward declaration
template<typename MeshType, typename BasisFunctionType>
class BasisOnMesh;
 
/** This adds functionality to create a partition / domain decomposition
 */
template<typename MeshType,typename BasisFunctionType>
class BasisOnMeshPartitionBase :
  public BasisOnMeshJacobian<MeshType,BasisFunctionType>
{
public:
  
  //! constructor
  BasisOnMeshPartitionBase(std::shared_ptr<Partition::Manager> partitionManager, PyObject *specificSettings);
  
  //! get the partition
  std::shared_ptr<Partition::MeshPartition<BasisOnMesh<MeshType,BasisFunctionType>,MeshType>> meshPartition() const;
  
  //! get the partition as pointer of type meshPartitionBase, this is in the itnerface in mesh
  std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBase();

protected:
  std::shared_ptr<Partition::Manager> partitionManager_;  ///< the partition manager object that can create partitions
  std::shared_ptr<Partition::MeshPartition<BasisOnMesh<MeshType,BasisFunctionType>,MeshType>> meshPartition_;   ///< the partition information that is stored locally, i.e. the subdomain of the domain decomposition

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
class BasisOnMeshPartition<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType> :
  public BasisOnMeshPartitionBase<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>
{
public:
  //! use inherited constructor
  using BasisOnMeshPartitionBase<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::BasisOnMeshPartitionBase;
  
  //! initiate the partitoning and then call the downwards initialize
  void initialize();
  
  //! get the total number of elements on the global domain, for structured meshes this is directly implemented in the Mesh itself (not BasisOnMesh like here)
  virtual global_no_t nElementsGlobal() const = 0;

};

}  // namespace

#include "basis_on_mesh/03_basis_on_mesh_partition_base.tpp"
#include "basis_on_mesh/03_basis_on_mesh_partition_structured.tpp"
#include "basis_on_mesh/03_basis_on_mesh_partition_unstructured.tpp"
