#pragma once

#include <memory>

#include "control/types.h"
#include "partition/rank_subset.h"
#include "partition/partitioned_petsc_mat_base.h"
#include "partition/01_mesh_partition.h"

/** Global numbering: such that each rank has its own contiguous subset of the total range.
 *  Local numbering: ghost elements
 */
template<typename BasisOnMeshType, typename = typename BasisOnMeshType::Mesh>
class PartitionedPetscMat
{
};

/** partial specialization for structured meshes */
template<typename MeshType, typename BasisFunctionType>
class PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>> : 
  public PartitionedPetscMatBase<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>>
{
public:
  //! constructor
  PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>>> meshPartition, 
                      int nComponents, int diagonalNonZeros, int offdiagonalNonZeros, std::string name);
 
protected:
  
  //! create a distributed Petsc matrix, according to the given partition
  void createMatrix(int diagonalNonZeros, int offdiagonalNonZeros);
  
  int nComponents_;  ///< number of components of the field variable
};

/** partial specialization for unstructured meshes 
 *  use petsc IS
 *  not implemented yet
 */
template<int D, typename BasisFunctionType>
class PartitionedPetscMat<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>> : 
  public PartitionedPetscMatBase<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>
{
public:
  //! constructor
  PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>> meshPartition, 
                      int nComponents, int diagonalNonZeros, int offdiagonalNonZeros, std::string name);
 
protected:
  
  //! create a distributed Petsc matrix, according to the given partition
  void createMatrix(int diagonalNonZeros, int offdiagonalNonZeros);
  
  int nComponents_;  ///< number of components of the field variable
  
};

template<typename BasisOnMeshType>
std::ostream &operator<<(std::ostream &stream, const PartitionedPetscMat<BasisOnMeshType> &matrix);

#include "partition/partitioned_petsc_mat_structured.tpp"
#include "partition/partitioned_petsc_mat_unstructured.tpp"
