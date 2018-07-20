#pragma once

#include <memory>

#include "control/types.h"
#include "partition/rank_subset.h"

/** This encapsulates a Petsc Vec, combined with the partition of the mesh.
 *  A local Vec is stored that holds all values 
 *  Global numbering: such that each rank has its own contiguous subset of the total range.
 *  Local numbering: including ghost elements
 */
template<typename BasisOnMeshType, int nComponents, typename = typename BasisOnMeshType::Mesh>
class PartitionedPetscVec
{
};

/** partial specialization for structured meshes */
template<int D, typename MeshType, typename BasisFunctionType, int nComponents>
class PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,nComponents,Mesh::isStructuredWithDim<D,MeshType>>
{
public:
 
  //! constructor, construct a petsc Vec with meshPartition that can hold the values for a field variable with nComponents components
  PartitionedPetscVec(MeshPartition &meshPartition, std::string name);
 
  //! this has to be called before the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called)
  void startVectorManipulation();
  
  //! this has to be called after the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called)
  void finishVectorManipulation();
  
  //! wrapper to the PETSc VecSetValues, acting only on the local data, the indices ix are the local dof nos
  void setValues(int componentNo, PetscInt ni, const PetscInt ix[], const PetscScalar y[], InsertMode iora);
  
  //! wrapper to the PETSc VecSetValue, acting only on the local data
  void setValue(int componentNo, PetscInt row, PetscScalar value, InsertMode mode);
  
  //! wrapper to the PETSc VecGetValues, acting only on the local data, the indices ix are the local dof nos
  void getValues(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[]);
  
  //! wrapper to the PETSc VecGetValues, on the global vector with global indexing
  void getValuesGlobalIndexing(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[]);
  
  //! set all entries to zero, wraps VecZeroEntries
  void zeroEntries();
  
  //! get the local Vector of a specified component
  Vec &values(int componentNo = 0);
  
protected:
 
  //! create a distributed Petsc vector, according to partition
  void createVector(std::string name);
  
  MeshPartition &meshPartition_;  ///< the mesh partition object which stores how the mesh is decomposed and what is the local portion
  DM dm_;    ///< PETSc DMDA object (nowhere specified what the abbreviation means) that stores topology information and everything needed for communication of ghost values
  
  std::array<Vec,nComponents> vectorLocal_;   ///< local vector that holds the local vectors, is filled by startVectorManipulation and can the be manipulated, afterwards the results need to get copied back by finishVectorManipulation
  std::array<Vec,nComponents> vectorGlobal_;  ///< the global distributed vector that holds the actual data
};

/** partial specialization for unstructured meshes 
 *  use petsc IS
 *  not implemented yet
 */
template<int D, typename BasisFunctionType, int nComponents>
class PartitionedPetscVec<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, nComponents, Mesh::UnstructuredDeformableOfDimension<D>>
{
public:
  PartitionedPetscVec(MeshPartition &meshPartition, std::string name, dof_no_t nEntries);
  
};

#include "partition/partitioned_petsc_vec.tpp"
