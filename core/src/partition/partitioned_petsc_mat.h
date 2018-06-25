#pragma once

#include <memory>

#include "control/types.h"
#include "partition/rank_subset.h"

/** Global numbering: such that each rank has its own contiguous subset of the total range.
 *  Local numbering: ghost elements
 */
template<typename BasisOnMeshType, typename = typename BasisOnMeshType::Mesh>
class PartitionedPetscMat
{
};

/** partial specialization for structured meshes */
template<int D, typename MeshType, typename BasisFunctionType>
class PartitionedPetscMat<BasisOnMesh::BasisOnMesh<MeshType,BasisFunctionType>,Mesh::isStructuredWithDim<D,MeshType>>
{
public:
 
  //! constructor
  PartitionedPetscMat(std::shared_ptr<MeshPartition> meshPartition, int nComponents, int diagonalNonZeros, int offdiagonalNonZeros);
 
  //! constructor to simply wrap an existing Mat, as needed in nonlinear solver callback functions for jacobians
  PartitionedPetscMat(Mat &matrix);
 
  //! wrapper of MatSetValues for a single value, sets a local value in the matrix
  void setValue(PetscInt row, PetscInt col, PetscScalar value, InsertMode mode);
  
  //! wrapper of MatZeroRowsColumns, zeros all entries (except possibly the main diagonal) of a set of local rows and columns
  void zeroRowsColumns(PetscInt numRows,const PetscInt rows[], PetscScalar diag);
  
  //! parallel assembly of the matrix, wraps the PETSc function MatAssemblyBegin,MatAssemblyEnd, type is MAT_FLUSH_ASSEMBLY or MAT_FINAL_ASSEMBLY
  void assembly(MatAssemblyType type);

  //! get entries from the matrix that are locally stored, uses the global indexing  
  void getValuesGlobalIndexing(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[]);
  
  //! get a reference to the PETSc matrix
  Mat &values();
  
protected:
 
  //! create a distributed Petsc matrix, according to the given partition
  void createMatrix();
  
  std::shared_ptr<MeshPartition> meshPartition_;  ///< the mesh partition object which stores how the mesh is decomposed and what is the local portion
  int nComponents_;  ///< number of components of the field variable
  
  Mat matrix_;   ///< the global Petsc matrix, access it performed using MatSetValuesLocal() with local indices
};

/** partial specialization for unstructured meshes 
 *  use petsc IS
 *  not implemented yet
 */
template<int D, typename BasisFunctionType>
class PartitionedPetscMat<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>>
{
public:
  
};

#include "partition/partitioned_petsc_vec.tpp"
