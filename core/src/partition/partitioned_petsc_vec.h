#pragma once

#include <memory>

#include "control/types.h"
#include "partition/rank_subset.h"

/** This encapsulates a Petsc Vec, combined with the partition of the mesh.
 *  A local Vec is stored that holds all values 
 *  Global numbering: such that each rank has its own contiguous subset of the total range.
 *  Local numbering: including ghost elements
 */
template<typename BasisOnMeshType, typename = typename BasisOnMeshType::Mesh>
class PartitionedPetscVec
{
};

/** partial specialization for structured meshes */
template<int D, typename MeshType, typename BasisFunctionType>
class PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>
{
public:
 
  //! constructor, construct a petsc Vec with meshPartition that can hold the values for a field variable with nComponents components
  PartitionedPetscVec(MeshPartition &meshPartition, int nComponents);
 
  //! this has to be called before the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called)
  void startVectorManipulation();
  
  //! this has to be called after the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called)
  void finishVectorManipulation();
  
  //! wrapper to the PETSc VecSetValues, acting only on the local data
  void setValues(PetscInt ni, const PetscInt ix[], const PetscScalar y[], InsertMode iora);
  
  //! wrapper to the PETSc VecGetValues, acting only on the local data
  void getValues(PetscInt ni, const PetscInt ix[], PetscScalar y[]);
  
  //! wrapper to the PETSc VecGetValues, on the global vector with global indexing
  void getValuesGlobalIndexing(PetscInt ni, const PetscInt ix[], PetscScalar y[]);
  
  //! set all entries to zero, wraps VecZeroEntries
  void zeroEntries();
  
  //! get the local Vector
  Vec &values();
  
protected:
 
  //! create a distributed Petsc vector, according to partition
  void createVector(std::string name);
  
  MeshPartition &meshPartition_;  ///< the mesh partition object which stores how the mesh is decomposed and what is the local portion
  DM dm_;    ///< PETSc DMDA object (nowhere specified what the abbreviation means) that stores topology information and everything needed for communication of ghost values
  int nComponents_;  ///< number of components of the field variable
  
  Vec vectorLocal_;   ///< local vector that holds the local vectors, is filled by startVectorManipulation and can the be manipulated, afterwards the results need to get copied back by finishVectorManipulation
  Vec vectorGlobal_;  ///< the global distributed vector that holds the actual data
};

/** partial specialization for unstructured meshes 
 *  use petsc IS
 *  not implemented yet
 */
template<int D, typename BasisFunctionType>
class PartitionedPetscVec<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>
{
public:
  
};

#include "partition/partitioned_petsc_vec.tpp"
