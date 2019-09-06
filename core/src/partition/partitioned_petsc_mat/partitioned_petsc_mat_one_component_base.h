#pragma once

#include <Python.h>  // has to be the first included header
#include <memory>

#include "control/types.h"
#include "partition/rank_subset.h"
#include "partition/mesh_partition/01_mesh_partition.h"

/** Base class for a partitioned PetscMat
 */
template<typename RowsFunctionSpaceType, typename ColumnsFunctionSpaceType>
class PartitionedPetscMatOneComponentBase
{
public:
  
  typedef RowsFunctionSpaceType RowsFunctionSpace;
  typedef ColumnsFunctionSpaceType ColumnsFunctionSpace;
  
  //! constructor
  PartitionedPetscMatOneComponentBase(std::shared_ptr<Partition::MeshPartition<RowsFunctionSpaceType>> meshPartitionRows,
                          std::shared_ptr<Partition::MeshPartition<ColumnsFunctionSpaceType>> meshPartitionColumns, std::string name);
 
  //! wrapper of MatSetValues for a single value, sets a local value in the matrix
  virtual void setValue(PetscInt row, PetscInt col, PetscScalar value, InsertMode mode) = 0;
  
  //! wrapper of MatSetValues for a single value, sets a local value in the matrix
  virtual void setValues(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], const PetscScalar v[], InsertMode addv) = 0;
  
  //! wrapper of MatZeroRowsColumns, zeros all entries (except possibly the main diagonal) of a set of local rows and columns
  virtual void zeroRowsColumns(PetscInt numRows,const PetscInt rows[], PetscScalar diag) = 0;
  
  //! wrapper of MatZeroEntries, sets all entries to 0
  virtual void zeroEntries() = 0;
  
  //! parallel assembly of the matrix, wraps the PETSc function MatAssemblyBegin,MatAssemblyEnd, type is MAT_FLUSH_ASSEMBLY or MAT_FINAL_ASSEMBLY
  virtual void assembly(MatAssemblyType type) = 0;

  //! get entries from the matrix that are locally stored, uses the local indexing
  virtual void getValues(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[]) const = 0;

  //! get entries from the matrix that are locally stored, uses the global indexing
  virtual void getValuesGlobalPetscIndexing(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[]) const = 0;
  
  //! get a reference to the local PETSc matrix
  virtual Mat &valuesLocal() = 0;
  
  //! get a reference to the global PETSc matrix
  virtual Mat &valuesGlobal() = 0;
    
  //! output matrix to stream, the operator<< is also overloaded to use this method
  virtual void output(std::ostream &stream) const = 0;
  
  //! get the mesh partition of rows
  std::shared_ptr<Partition::MeshPartition<RowsFunctionSpaceType>> meshPartitionRows();
  
  //! get the mesh partion of columns
  std::shared_ptr<Partition::MeshPartition<ColumnsFunctionSpaceType>> meshPartitionColumns();
  
protected:
 
  std::shared_ptr<Partition::MeshPartition<RowsFunctionSpaceType>> meshPartitionRows_;  ///< the mesh partition object which stores how the mesh is decomposed and what is the local portion, for the rows of the matrix
  std::shared_ptr<Partition::MeshPartition<ColumnsFunctionSpaceType>> meshPartitionColumns_;  ///< the mesh partition object which stores how the mesh is decomposed and what is the local portion, for the columns of the matrix
  std::string name_;   ///< a specifier for the matrix, only used for debugging
};


#include "partition/partitioned_petsc_mat/partitioned_petsc_mat_one_component_base.tpp"
