#pragma once

#include <memory>
#include "partition/mesh_partition/01_mesh_partition_structured.h"
#include "partition/mesh_partition/01_mesh_partition_unstructured.h"
#include "partition/mesh_partition/01_mesh_partition_composite.h"
#include "partition/partitioned_petsc_vec/values_representation.h"

/** Base class for partitioned petsc vectors, this just holds a pointer to the meshPartion object and stores the data representation.
 * There are 3 representations:
 * representationLocal: This is the normal case, getValues and setValues work on the local vectors. Ghost dofs are included in the local vectors.
 * representationGlobal: The actual information is in the global vectors. This is needed, when the values should be accessed using valuesGlobal(), for direct use with Petsc functions.
 * representationContiguous: The current local data is stored in a contiguous way, i.e. all components after each other. This is needed for access by CellML.
 *
 */
template<typename FunctionSpaceType>
class PartitionedPetscVecBase
{
public:
 
  //! constructor
  PartitionedPetscVecBase(std::shared_ptr<Partition::MeshPartition<FunctionSpaceType,typename FunctionSpaceType::Mesh>> meshPartition, std::string name);
 
  //! get a vector of local dof nos (from meshPartition), with ghost dofs
  std::vector<PetscInt> &localDofNos();
  
  //! get the meshPartition
  const std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> meshPartition() const;
  
  //! get the name of the vector
  std::string name() const;

  //! get the current internal data representation
  Partition::values_representation_t currentRepresentation() const;

  //! get a string of the current representation
  std::string getCurrentRepresentationString() const;

protected:
  
  std::string name_;   //< name of the vector
  std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> meshPartition_;  //< the mesh partition object which stores how the mesh is decomposed and what is the local portion

  Partition::values_representation_t currentRepresentation_;  //< which vector holds the current values, valuesLocal, valuesGlobal or valuesContiguous

  static int vectorNo_;   //< a contiguous number of the vector which is added to the vector name
};


#include "partition/partitioned_petsc_vec/partitioned_petsc_vec_base.tpp"
