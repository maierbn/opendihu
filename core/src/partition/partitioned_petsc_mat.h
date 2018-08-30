#pragma once

#include <memory>

#include "control/types.h"
#include "partition/rank_subset.h"
#include "partition/partitioned_petsc_mat_base.h"
#include "partition/01_mesh_partition.h"

/** Global numbering: such that each rank has its own contiguous subset of the total range.
 *  Local numbering: ghost elements
 */
template<typename FunctionSpaceType, typename = typename FunctionSpaceType::Mesh>
class PartitionedPetscMat
{
};

/** partial specialization for structured meshes */
template<typename MeshType, typename BasisFunctionType>
class PartitionedPetscMat<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,Mesh::isStructured<MeshType>> : 
  public PartitionedPetscMatBase<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>
{
public:
  //! constructor
  PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>> meshPartition, 
                      int nComponents, int diagonalNonZeros, int offdiagonalNonZeros, std::string name);
 
  //! constructor, use provided global matrix
  PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>> meshPartition,
                      Mat &globalMatrix, std::string name);

  //! wrapper of MatSetValues for a single value, sets a local value in the matrix
  void setValue(PetscInt row, PetscInt col, PetscScalar value, InsertMode mode);
  
  //! wrapper of MatSetValues for a single value, sets a local value in the matrix
  void setValues(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], const PetscScalar v[], InsertMode addv);
  
  //! wrapper of MatZeroRowsColumns, zeros all entries (except possibly the main diagonal) of a set of local rows and columns
  void zeroRowsColumns(PetscInt numRows, const PetscInt rows[], PetscScalar diag);
  
  //! wrapper of MatZeroEntries, sets all entries to 0
  void zeroEntries();
  
  //! parallel assembly of the matrix, wraps the PETSc function MatAssemblyBegin,MatAssemblyEnd, type is MAT_FLUSH_ASSEMBLY or MAT_FINAL_ASSEMBLY
  void assembly(MatAssemblyType type);

  //! get entries from the matrix that are locally stored, uses the local indexing
  void getValues(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[]) const;

  //! get entries from the matrix that are locally stored, uses the global indexing
  void getValuesGlobalIndexing(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[]);

  //! get a reference to the local PETSc matrix
  Mat &valuesLocal();
  
  //! get a reference to the global PETSc matrix
  Mat &valuesGlobal();
    
  //! output matrix to stream, the operator<< is also overloaded to use this method
  void output(std::ostream &stream) const;
  
protected:
  
  //! create a distributed Petsc matrix, according to the given partition
  void createMatrix(int diagonalNonZeros, int offdiagonalNonZeros);

  //! set the global to local mapping at the global matrix and create the local submatrix
  void createLocalMatrix();

  Mat globalMatrix_;   ///< the global Petsc matrix, access using MatSetValuesLocal() with local indices (not used here) or via the localMatrix (this one is used)
  Mat localMatrix_;    ///< a local submatrix that holds all rows and columns for the local dofs with ghosts
  
  int nComponents_;  ///< number of components of the field variable
};

/** Partial specialization for unstructured meshes. This is completely serial, there are no parallel matrices.
 *  (To enable for parallelism, use petsc IS, not implemented yet)
 */
template<int D, typename BasisFunctionType>
class PartitionedPetscMat<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>, Mesh::UnstructuredDeformableOfDimension<D>> : 
  public PartitionedPetscMatBase<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>
{
public:
  //! constructor, create matrix
  PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>> meshPartition,
                      int nComponents, int diagonalNonZeros, int offdiagonalNonZeros, std::string name);

  //! constructor, use provided global matrix
  PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>> meshPartition,
                      Mat &globalMatrix, std::string name);

  //! wrapper of MatSetValues for a single value, sets a local value in the matrix
  void setValue(PetscInt row, PetscInt col, PetscScalar value, InsertMode mode);
  
  //! wrapper of MatSetValues for a single value, sets a local value in the matrix
  void setValues(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], const PetscScalar v[], InsertMode addv);
  
  //! wrapper of MatZeroRowsColumns, zeros all entries (except possibly the main diagonal) of a set of local rows and columns
  void zeroRowsColumns(PetscInt numRows,const PetscInt rows[], PetscScalar diag);
  
  //! wrapper of MatZeroEntries, sets all entries to 0
  void zeroEntries();
  
  //! parallel assembly of the matrix, wraps the PETSc function MatAssemblyBegin,MatAssemblyEnd, type is MAT_FLUSH_ASSEMBLY or MAT_FINAL_ASSEMBLY
  void assembly(MatAssemblyType type);

  //! get entries from the matrix that are locally stored
  void getValues(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[]) const;

  //! get entries from the matrix that are locally stored
  void getValuesGlobalIndexing(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[]);
  
  //! get a reference to the PETSc matrix, because there is no parallelism with UnstructuredDeformableOfDimension meshes, this is the same as valuesGlobal
  Mat &valuesLocal();
  
  //! get a reference to the PETSc matrix, because there is no parallelism with UnstructuredDeformableOfDimension meshes, this is the same as valuesLocal
  Mat &valuesGlobal();
    
  //! output matrix to stream, the operator<< is also overloaded to use this method
  void output(std::ostream &stream) const;
  
protected:
  
  //! create a distributed Petsc matrix, according to the given partition
  void createMatrix(int diagonalNonZeros, int offdiagonalNonZeros);
  
  Mat matrix_;   ///< the single Petsc matrix (global = local)
  int nComponents_;  ///< number of components of the field variable
  
};

template<typename FunctionSpaceType>
std::ostream &operator<<(std::ostream &stream, const PartitionedPetscMat<FunctionSpaceType> &matrix);

#include "partition/partitioned_petsc_mat_structured.tpp"
#include "partition/partitioned_petsc_mat_unstructured.tpp"
