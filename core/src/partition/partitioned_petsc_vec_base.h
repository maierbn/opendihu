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
  PartitionedPetscVecBase(std::shared_ptr<Partition::MeshPartition<BasisOnMeshType,typename BasisOnMeshType::Mesh>> meshPartition);
 
  //! get a vector of local dofs (from meshPartition)
  std::vector<PetscInt> &localDofs();
  
protected:
  
  std::shared_ptr<Partition::MeshPartition<BasisOnMeshType>> meshPartition_;  ///< the mesh partition object which stores how the mesh is decomposed and what is the local portion
};

#include "partition/partitioned_petsc_vec_base.tpp"
