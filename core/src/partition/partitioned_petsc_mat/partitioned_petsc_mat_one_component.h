#pragma once

#include <Python.h>  // has to be the first included header
#include <memory>

#include "control/types.h"
#include "partition/rank_subset.h"
//#include "partition/partitioned_petsc_mat/partitioned_petsc_mat.h"
#include "partition/partitioned_petsc_mat/partitioned_petsc_mat_one_component_base.h"
#include "partition/mesh_partition/01_mesh_partition.h"

/** Global numbering: such that each rank has its own contiguous subset of the total range.
 *  Local numbering: ghost elements
 */
template<typename RowsFunctionSpaceType, typename ColumnsFunctionSpaceType = RowsFunctionSpaceType>
class PartitionedPetscMatOneComponent
{
};

/** partial specialization for structured meshes and composite mesh*/
template<typename MeshType, typename BasisFunctionType, typename ColumnsFunctionSpaceType>
class PartitionedPetscMatOneComponent<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,ColumnsFunctionSpaceType> :
  public PartitionedPetscMatOneComponentBase<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>,ColumnsFunctionSpaceType>
{
public:
  //! constructor, create square sparse matrix
  PartitionedPetscMatOneComponent(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>> meshPartition,
                                  int nNonZerosDiagonal, int nNonZerosOffdiagonal, std::string name);

  //! constructor, create square dense matrix
  PartitionedPetscMatOneComponent(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>> meshPartition,
                                  std::string name);

  //! constructor, create non-square sparse matrix
  PartitionedPetscMatOneComponent(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>> meshPartitionRows,
                                  std::shared_ptr<Partition::MeshPartition<ColumnsFunctionSpaceType>> meshPartitionColumns,
                                  int nNonZerosDiagonal, int nNonZerosOffdiagonal, std::string name);

  //! constructor, create non-square dense matrix
  PartitionedPetscMatOneComponent(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>> meshPartitionRows,
                                  std::shared_ptr<Partition::MeshPartition<ColumnsFunctionSpaceType>> meshPartitionColumns,
                                  std::string name);

  //! constructor, use provided global matrix
  PartitionedPetscMatOneComponent(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<MeshType,BasisFunctionType>>> meshPartition,
                                  Mat &globalMatrix, std::string name);

  //! destructor
  virtual ~PartitionedPetscMatOneComponent();

  //! wrapper of MatSetValues for a single value, sets a local value in the matrix
  void setValue(PetscInt row, PetscInt col, PetscScalar value, InsertMode mode);
  
  //! wrapper of MatSetValues for a single value, sets a local value in the matrix
  void setValues(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], const PetscScalar v[], InsertMode addv);
  
  //! set entries in the matrix that are locally stored, uses the global/Petsc indexing. This is not the global natural numbering!
  void setValuesGlobalPetscIndexing(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], const PetscScalar v[], InsertMode addv);

  //! wrapper of MatZeroRowsColumns, zeros all entries (except possibly the main diagonal) of a set of local rows and columns
  void zeroRowsColumns(PetscInt numRows, const PetscInt rows[], PetscScalar diag);

  //! wrapper of MatZeroRows, zeros all entries local rows
  void zeroRows(PetscInt numRows, const PetscInt rows[]);

  //! wrapper of MatZeroEntries, sets all entries to 0
  void zeroEntries();
  
  //! parallel assembly of the matrix, wraps the PETSc function MatAssemblyBegin,MatAssemblyEnd, type is MAT_FLUSH_ASSEMBLY or MAT_FINAL_ASSEMBLY
  void assembly(MatAssemblyType type);

  //! get entries from the matrix that are locally stored, uses the local indexing
  void getValues(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[]) const;

  //! get entries from the matrix that are locally stored, uses the global/Petsc indexing. This is not the global natural numbering!
  void getValuesGlobalPetscIndexing(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[]) const;

  //! get a reference to the local PETSc matrix
  Mat &valuesLocal();
  
  //! get a reference to the global PETSc matrix
  Mat &valuesGlobal();
    
  //! output matrix to stream, the operator<< is also overloaded to use this method
  void output(std::ostream &stream) const;
  
  //! write the vector to a file using PetscViewer, format is "default", "ascii" or "matlab"
  void dumpMatrix(std::string filename, std::string format);

protected:
  
  //! create a distributed Petsc matrix, according to the given partition
  void createMatrix(MatType matrixType, int nNonZerosDiagonal, int nNonZerosOffdiagonal);

  //! set the global to local mapping at the global matrix and create the local submatrix
  void createLocalMatrix();

  Mat globalMatrix_;   ///< the global Petsc matrix, access using MatSetValuesLocal() with local indices (not used here) or via the localMatrix (this one is used)
  Mat localMatrix_;    ///< a local submatrix that holds all rows and columns for the local dofs with ghosts
};

/** Partial specialization for unstructured meshes. This is completely serial, there are no parallel matrices.
 *  (To enable for parallelism, use petsc IS, not implemented yet)
 */
template<int D, typename BasisFunctionType>
class PartitionedPetscMatOneComponent<
  FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>,
  FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>> :
  public PartitionedPetscMatOneComponentBase<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>,FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunctionType>>
{
public:
  //! constructor, create square sparse matrix
  PartitionedPetscMatOneComponent(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>> meshPartition,
                                  int nNonZerosDiagonal, int nNonZerosOffdiagonal, std::string name);

  //! constructor, create square dense matrix
  PartitionedPetscMatOneComponent(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>> meshPartition,
                                  std::string name);

  //! constructor, create non-square sparse matrix
  PartitionedPetscMatOneComponent(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>> meshPartitionRows,
                                  std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>> meshPartitionColumns,
                                  int nNonZerosDiagonal, int nNonZerosOffdiagonal, std::string name);

  //! constructor, create non-square dense matrix
  PartitionedPetscMatOneComponent(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>> meshPartitionRows,
                                  std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>> meshPartitionColumns,
                                  std::string name);

  //! constructor, use provided global matrix
  PartitionedPetscMatOneComponent(std::shared_ptr<Partition::MeshPartition<FunctionSpace::FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>> meshPartition,
                                  Mat &globalMatrix, std::string name);

  //! wrapper of MatSetValues for a single value, sets a local value in the matrix
  void setValue(PetscInt row, PetscInt col, PetscScalar value, InsertMode mode);
  
  //! wrapper of MatSetValues for a single value, sets a local value in the matrix
  void setValues(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], const PetscScalar v[], InsertMode addv);

  //! set entries in the matrix that are locally stored, uses the global/Petsc indexing. This is not the global natural numbering!
  void setValuesGlobalPetscIndexing(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], const PetscScalar v[], InsertMode addv);

  //! wrapper of MatZeroRowsColumns, zeros all entries (except possibly the main diagonal) of a set of local rows and columns
  void zeroRowsColumns(PetscInt numRows,const PetscInt rows[], PetscScalar diag);
  
  //! wrapper of MatZeroRows, zeros all entries local rows
  void zeroRows(PetscInt numRows, const PetscInt rows[]);

  //! wrapper of MatZeroEntries, sets all entries to 0
  void zeroEntries();
  
  //! parallel assembly of the matrix, wraps the PETSc function MatAssemblyBegin,MatAssemblyEnd, type is MAT_FLUSH_ASSEMBLY or MAT_FINAL_ASSEMBLY
  void assembly(MatAssemblyType type);

  //! get entries from the matrix that are locally stored
  void getValues(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[]) const;

  //! get entries from the matrix that are locally stored, uses the global/Petsc indexing. This is not the global natural numbering!
  void getValuesGlobalPetscIndexing(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[]) const;
  
  //! get a reference to the PETSc matrix, because there is no parallelism with UnstructuredDeformableOfDimension meshes, this is the same as valuesGlobal
  Mat &valuesLocal();
  
  //! get a reference to the PETSc matrix, because there is no parallelism with UnstructuredDeformableOfDimension meshes, this is the same as valuesLocal
  Mat &valuesGlobal();
    
  //! output matrix to stream, the operator<< is also overloaded to use this method
  void output(std::ostream &stream) const;

  //! write the vector to a file using PetscViewer, format is "default", "ascii" or "matlab"
  void dumpMatrix(std::string filename, std::string format);

protected:
  
  //! create a distributed Petsc matrix, according to the given partition
  void createMatrix(MatType matrixType, int nNonZerosDiagonal, int nNonZerosOffdiagonal);
  
  Mat matrix_;   ///< the single Petsc matrix (global = local)
};

template<typename FunctionSpaceType>
std::ostream &operator<<(std::ostream &stream, const PartitionedPetscMatOneComponent<FunctionSpaceType> &matrix);

#include "partition/partitioned_petsc_mat/partitioned_petsc_mat_one_component_structured.tpp"
#include "partition/partitioned_petsc_mat/partitioned_petsc_mat_one_component_unstructured.tpp"
