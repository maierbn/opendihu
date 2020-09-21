#pragma once

#include <memory>

#include "control/types.h"
#include "partition/rank_subset.h"
#include "partition/mesh_partition/01_mesh_partition.h"
#include "partition/partitioned_petsc_mat/partitioned_petsc_mat_one_component.h"

/** Class that represents a Petsc Mat. This is potentially a nested matrix with multiple submatrices, if nComponents != 1.
 */
template<typename RowsFunctionSpaceType, typename ColumnsFunctionSpaceType = RowsFunctionSpaceType>
class PartitionedPetscMat
{
public:

  //! constructor, create square sparse matrix
  PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<RowsFunctionSpaceType>> meshPartition,
                      int nComponents, int nNonZerosDiagonal, int nNonZerosOffdiagonal, std::string name);

  //! constructor, create square dense matrix
  PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<RowsFunctionSpaceType>> meshPartition,
                      int nComponents, std::string name);

  //! constructor, create non-square sparse matrix
  PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<RowsFunctionSpaceType>> meshPartitionRows,
                      std::shared_ptr<Partition::MeshPartition<ColumnsFunctionSpaceType>> meshPartitionColumns,
                      int nComponents, int nNonZerosDiagonal, int nNonZerosOffdiagonal, std::string name);

  //! constructor, create non-square dense matrix
  PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<RowsFunctionSpaceType>> meshPartitionRows,
                      std::shared_ptr<Partition::MeshPartition<ColumnsFunctionSpaceType>> meshPartitionColumns,
                      int nComponents, std::string name);

  //! constructor, use provided global matrix
  PartitionedPetscMat(std::shared_ptr<Partition::MeshPartition<RowsFunctionSpaceType>> meshPartition,
                      Mat &globalMatrix, std::string name);

  //! wrapper of MatSetValues for a single value or component 0, sets a local value in the matrix
  void setValue(PetscInt row, PetscInt col, PetscScalar value, InsertMode mode);

  //! wrapper of MatSetValues for a single value, sets a local value in the matrix
  void setValue(int componentNo, PetscInt row, PetscInt col, PetscScalar value, InsertMode mode);

  //! wrapper of MatSetValues for a single value, sets a local value in the matrix
  void setValue(int componentNo, Vc::int_v row, Vc::int_v columns, PetscScalar value, InsertMode mode);

  //! wrapper of MatSetValues for a single value, sets a local value in the matrix
  void setValue(int componentNo, Vc::int_v row, Vc::int_v columns, Vc::double_v value, InsertMode mode);

  //! wrapper of MatSetValues for a single value, sets a local value in the matrix
  template<int nComponents>
  void setValue(PetscInt row, PetscInt col, std::array<double,nComponents> value, InsertMode mode);

  //! wrapper of MatSetValues for a single value or component 0, sets a local value in the matrix
  void setValues(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], const PetscScalar v[], InsertMode addv);

  //! wrapper of MatSetValues for a single value, sets a local value in the matrix
  void setValues(int componentNo, PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], const PetscScalar v[], InsertMode addv);

  //! wrapper of MatSetValues for a single value, sets a local value in the matrix
  template<int nComponents>
  void setValues(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], const std::vector<std::array<double,nComponents>> &v, InsertMode addv);

  //! set entries in the given submatrix, uses the global/Petsc indexing. This is not the global natural numbering!
  void setValuesGlobalPetscIndexing(int componentNo, PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], const PetscScalar v[], InsertMode addv);

  //! wrapper of MatZeroRowsColumns, zeros all entries (except possibly the main diagonal) of a set of local rows and columns
  void zeroRowsColumns(PetscInt numRows, const PetscInt rows[], PetscScalar diag);

  //! wrapper of MatZeroRowsColumns, zeros all entries (except possibly the main diagonal) of a set of local rows and columns
  //! this is for a given row / column component
  void zeroRowsColumns(int rowColumncomponentNo, PetscInt numRows, const PetscInt rows[], PetscScalar diag);

  //! wrapper of MatZeroEntries, sets all entries to 0
  void zeroEntries();

  //! parallel assembly of the matrix, wraps the PETSc function MatAssemblyBegin,MatAssemblyEnd, type is MAT_FLUSH_ASSEMBLY or MAT_FINAL_ASSEMBLY
  void assembly(MatAssemblyType type);

  //! get entries from the matrix that are locally stored
  void getValues(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[]) const;

  //! get entries from the matrix that are locally stored
  void getValues(int componentNo, PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[]) const;

  //! get entries from the matrix that are locally stored
  template<int nComponents>
  void getValues(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], std::vector<std::array<double,nComponents>> &v) const;

  //! get entries from the matrix that are locally stored, uses the global/Petsc indexing. This is not the global natural numbering!
  void getValuesGlobalPetscIndexing(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[]) const;

  //! get entries from the matrix that are locally stored, uses the global/Petsc indexing. This is not the global natural numbering!
  void getValuesGlobalPetscIndexing(int componentNo, PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[]) const;

  //! get entries from the matrix that are locally stored, uses the global/Petsc indexing. This is not the global natural numbering!
  template<int nComponents>
  void getValuesGlobalPetscIndexing(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], std::vector<std::array<double,nComponents>> &v) const;

  //! get a reference to the local PETSc matrix
  Mat &valuesLocal(int componentNo = 0);

  //! get a reference to the global PETSc matrix of the specified component
  Mat &valuesGlobal(int componentNo);

  //! get a reference to the nested global PETSc matrix
  Mat &valuesGlobal();

  //! output matrix to stream, the operator<< is also overloaded to use this method
  void output(std::ostream &stream) const;

  //! get the mesh partition of rows
  std::shared_ptr<Partition::MeshPartition<RowsFunctionSpaceType>> meshPartitionRows();

  //! get the mesh partion of columns
  std::shared_ptr<Partition::MeshPartition<ColumnsFunctionSpaceType>> meshPartitionColumns();

  //! write the matrix to a file using PetscViewer, format is "default", "ascii" or "matlab"
  void dumpMatrix(std::string filename, std::string format);

protected:

  //! create the nested matrix
  void createMatNest();

  std::vector<PartitionedPetscMatOneComponent<RowsFunctionSpaceType,ColumnsFunctionSpaceType>> matrixComponents_; // submatrices for nComponents_^2 (row major)
  Mat matNest_;    //< nested matrix of MatNest type that contains the submatrices

  int nComponents_;  //< number of components of the field variable that can be right-multiplied with this matrix, number of submatrices is nComponents_^2
};

template<typename FunctionSpaceType>
std::ostream &operator<<(std::ostream &stream, const PartitionedPetscMat<FunctionSpaceType> &matrix);

#include "partition/partitioned_petsc_mat/partitioned_petsc_mat.tpp"
