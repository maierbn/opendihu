#pragma once

#include <Python.h>  // has to be the first included header

#include "function_space/02_function_space_jacobian.h"
#include "partition/01_mesh_partition.h"
#include "partition/00_mesh_partition_base.h"

// forward declaration
namespace Partition 
{
class Manager;
//template<typename FunctionSpaceType, typename DummyForTraits>
//class MeshPartition;
}

namespace FunctionSpace
{
 
// forward declaration
template<typename MeshType, typename BasisFunctionType>
class FunctionSpace;
 
/** This adds functionality to create a partition / domain decomposition
 */
template<typename MeshType,typename BasisFunctionType>
class FunctionSpacePartitionBase :
  public FunctionSpaceJacobian<MeshType,BasisFunctionType>
{
public:
  
  //! constructor
  FunctionSpacePartitionBase(std::shared_ptr<Partition::Manager> partitionManager, PyObject *specificSettings);
  
  //! get the partition
  std::shared_ptr<Partition::MeshPartition<FunctionSpace<MeshType,BasisFunctionType>,MeshType>> meshPartition() const;
  
  //! get the partition as pointer of type meshPartitionBase, this is in the itnerface in mesh
  std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBase();

protected:
  std::shared_ptr<Partition::Manager> partitionManager_;  ///< the partition manager object that can create partitions
  std::shared_ptr<Partition::MeshPartition<FunctionSpace<MeshType,BasisFunctionType>,MeshType>> meshPartition_;   ///< the partition information that is stored locally, i.e. the subdomain of the domain decomposition

};

/** specialization for structured meshes 
 */
template<typename MeshType,typename BasisFunctionType>
class FunctionSpacePartition :
  public FunctionSpacePartitionBase<MeshType,BasisFunctionType>
{
public:
  //! use inherited constructor 
  using FunctionSpacePartitionBase<MeshType,BasisFunctionType>::FunctionSpacePartitionBase;
  
  //! initiate the partitoning and then call the downwards initialize
  void initialize();

};

/** partial specialization for unstructured mesh
 */
template<int D,typename BasisFunctionType>
class FunctionSpacePartition<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType> :
  public FunctionSpacePartitionBase<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>
{
public:
  //! use inherited constructor
  using FunctionSpacePartitionBase<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::FunctionSpacePartitionBase;
  
  //! initiate the partitoning and then call the downwards initialize
  void initialize();
  
  //! get the total number of elements on the global domain, for structured meshes this is directly implemented in the Mesh itself (not FunctionSpace like here)
  virtual global_no_t nElementsGlobal() const = 0;

  //! get the number of nodes, this is stored in geometry field
  virtual global_no_t nNodesGlobal() const = 0;
  
  //! get the number of dofs, this is stored in this->nDofs_ and in geometry field
  virtual global_no_t nDofsGlobal() const = 0;

  // nDofsGlobal() is defined in 06_function_space_dofs_nodes.h
  
};

}  // namespace

#include "function_space/03_function_space_partition_base.tpp"
#include "function_space/03_function_space_partition_structured.tpp"
#include "function_space/03_function_space_partition_unstructured.tpp"
