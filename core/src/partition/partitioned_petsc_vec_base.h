#pragma once

#include <memory>
#include "partition/01_mesh_partition.h"

/** Base class for partitioned petsc vectors, this just holds a pointer to the meshPartion object.
 */
template<typename BasisOnMeshType>
class PartitionedPetscVecBase
{
public:
 
  //! constructor
  PartitionedPetscVecBase(std::shared_ptr<Partition::MeshPartition<BasisOnMeshType,typename BasisOnMeshType::Mesh>> meshPartition, std::string name);
 
  //! get a vector of local dof nos (from meshPartition), with ghost dofs
  std::vector<PetscInt> &localDofNos();
  
protected:
  
  std::string name_;   ///< name of the vector
  std::shared_ptr<Partition::MeshPartition<BasisOnMeshType>> meshPartition_;  ///< the mesh partition object which stores how the mesh is decomposed and what is the local portion
};

#include "partition/partitioned_petsc_vec_base.tpp"
