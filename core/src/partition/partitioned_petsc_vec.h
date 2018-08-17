#pragma once

#include <memory>

#include "control/types.h"
#include "partition/rank_subset.h"
#include "basis_on_mesh/basis_on_mesh.h"
#include "mesh/type_traits.h"
#include "partition/partitioned_petsc_vec_base.h"

// forward declaration
namespace BasisOnMesh
{
template<typename MeshType, typename BasisFunctionType>
class BasisOnMesh;
};

/** This encapsulates a Petsc Vec, combined with the partition of the mesh.
 *  A local Vec is stored that holds all values 
 *  Global numbering: such that each rank has its own contiguous subset of the total range.
 *  Local numbering: starting with 0, including ghost elements
 * *
 *  This particular standard specialization is for non-structured meshes or no meshes and currently completely serial, 
 *  it is the placeholder as long as the partial specialization for unstructured meshes is not implemented.
 */
template<typename BasisOnMeshType, int nComponents, typename = typename BasisOnMeshType::Mesh>
class PartitionedPetscVec : 
  public PartitionedPetscVecBase<BasisOnMeshType>
{
public:
  //! constructor
  PartitionedPetscVec(std::shared_ptr<Partition::MeshPartition<BasisOnMeshType,typename BasisOnMeshType::Mesh>> meshPartition, std::string name);
 
  //! constructor, copy from existing vector
  PartitionedPetscVec(PartitionedPetscVec<BasisOnMeshType,nComponents> &rhs, std::string name);
  
  //! this has to be called before the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called)
  void startVectorManipulation();
  
  //! this has to be called after the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called)
  void finishVectorManipulation();
  
  //! wrapper to the PETSc VecSetValues, acting only on the local data, the indices ix are the local dof nos
  void setValues(int componentNo, PetscInt ni, const PetscInt ix[], const PetscScalar y[], InsertMode iora);
  
  //! wrapper to the PETSc VecSetValue, acting only on the local data
  void setValue(int componentNo, PetscInt row, PetscScalar value, InsertMode mode);
  
  //! for a single component vector set all values, the input vector values is expected to have ghosts included
  void setValuesWithGhosts(int componentNo, std::vector<double> &values, InsertMode petscInsertMode);
  
  //! for a single component vector set all values, the input vector values is expected to have no ghosts included
  void setValuesWithoutGhosts(int componentNo, std::vector<double> &values, InsertMode petscInsertMode);
  
  //! wrapper to the PETSc VecGetValues, acting only on the local data, the indices ix are the local dof nos
  void getValues(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[]);
  
  //! wrapper to the PETSc VecGetValues, on the global vector with global indexing
  void getValuesGlobalIndexing(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[]);
  
  //! get all locally stored values
  void getLocalValues(int componentNo, std::vector<double> &values);
  
  //! set all entries to zero, wraps VecZeroEntries
  void zeroEntries();
  
  //! get the local Vector of a specified component
  Vec &values(int componentNo = 0);

protected:
  
  //! create the values vectors
  void createVector();
  
  std::array<Vec,nComponents> values_;  // the (serial) Petsc vector that contains all the data
};

/** This is the partial specialization for structured meshes.
 *  An own DMDA object is generated, separately from the one in MeshPartition. This object now refers to nodes (as opposite the one of MeshPartition which refers to elements). 
 *  This object is created such that it matches the partition given by the meshPartition.
 */
template<typename MeshType, typename BasisFunctionType, int nComponents>
class PartitionedPetscVec<
  BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,
  nComponents,
  Mesh::isStructured<MeshType>> : 
  public PartitionedPetscVecBase<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>>
{
public:
 
  //! constructor, construct a petsc Vec with meshPartition that can hold the values for a field variable with nComponents components
  PartitionedPetscVec(std::shared_ptr<Partition::MeshPartition<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,MeshType>> meshPartition, std::string name);
 
  //! constructor, copy from existing vector
  PartitionedPetscVec(PartitionedPetscVec<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,nComponents> &rhs, std::string name);
 
  //! this has to be called before the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called)
  void startVectorManipulation();
  
  //! this has to be called after the vector is manipulated (i.e. VecSetValues or vecZeroEntries is called)
  void finishVectorManipulation();
  
  //! wrapper to the PETSc VecSetValues, acting only on the local data, the indices ix are the local dof nos
  void setValues(int componentNo, PetscInt ni, const PetscInt ix[], const PetscScalar y[], InsertMode iora);
  
  //! wrapper to the PETSc VecSetValue, acting only on the local data
  void setValue(int componentNo, PetscInt row, PetscScalar value, InsertMode mode);

  //! for a single component vector set all values. They have to be enough for all local dof including ghosts.
  void setValuesWithGhosts(int componentNo, std::vector<double> &values, InsertMode petscInsertMode);
  
  //! for a single component vector set all values. values does not contain ghost dofs.
  void setValuesWithoutGhosts(int componentNo, std::vector<double> &values, InsertMode petscInsertMode);

  //! wrapper to the PETSc VecGetValues, acting only on the local data, the indices ix are the local dof nos
  void getValues(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[]);
  
  //! wrapper to the PETSc VecGetValues, on the global vector with global indexing
  void getValuesGlobalIndexing(int componentNo, PetscInt ni, const PetscInt ix[], PetscScalar y[]);
  
  //! get all locally stored values, i.e. with ghosts
  void getLocalValues(int componentNo, std::vector<double> &values);
  
  //! set all entries to zero, wraps VecZeroEntries
  void zeroEntries();
  
  //! get the local Vector of a specified component
  Vec &values(int componentNo = 0);
  
  //! get a vector of local dof nos (from meshPartition), without ghost dofs
  std::vector<PetscInt> &localDofNosWithoutGhosts();
  
protected:
 
  //! create a distributed Petsc vector, according to partition
  void createVector();
  
  std::shared_ptr<DM> dm_;    ///< PETSc DMDA object that stores topology information and everything needed for communication of ghost values
  
  std::array<Vec,nComponents> vectorLocal_;   ///< local vector that holds the local vectors, is filled by startVectorManipulation and can the be manipulated, afterwards the results need to get copied back by finishVectorManipulation
  std::array<Vec,nComponents> vectorGlobal_;  ///< the global distributed vector that holds the actual data
};

/** partial specialization for unstructured meshes 
 *  use petsc IS
 *  not implemented yet, do not remove
 */
/*
template<int D, typename BasisFunctionType, int nComponents>
class PartitionedPetscVec<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, nComponents, Mesh::UnstructuredDeformableOfDimension<D>>
{
public:
  PartitionedPetscVec(MeshPartition &meshPartition, std::string name, dof_no_t nEntries);
  
 
};
*/

#include "partition/partitioned_petsc_vec_default.tpp"
#include "partition/partitioned_petsc_vec_structured.tpp"
