#pragma once

#include <memory>

#include "control/types.h"
#include "partition/rank_subset.h"

/** Base class for a partitioned PetscMat
 */
template<typename BasisOnMeshType>
class PartitionedPetscMatBase
{
public:
  
  //! constructor
  PartitionedPetscMatBase(std::shared_ptr<Partition::MeshPartition<BasisOnMeshType>> meshPartition, std::string name);
 
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

  //! get entries from the matrix that are locally stored, uses the global indexing  
  void getValuesGlobalIndexing(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[]);
  
  //! get a reference to the local PETSc matrix
  Mat &valuesLocal();
  
  //! get a reference to the global PETSc matrix
  Mat &valuesGlobal();
    
  //! output matrix to stream, the operator<< is also overloaded to use this method
  void output(std::ostream &stream) const;
  
protected:
 
  std::shared_ptr<Partition::MeshPartition<BasisOnMeshType>> meshPartition_;  ///< the mesh partition object which stores how the mesh is decomposed and what is the local portion
  Mat globalMatrix_;   ///< the global Petsc matrix, access using MatSetValuesLocal() with local indices (not used here) or via the localMatrix (this one is used)
  Mat localMatrix_;    ///< a local submatrix that holds all rows and columns for the local dofs with ghosts
  std::string name_;   ///< a specifier for the matrix, only used for debugging
};


#include "partition/partitioned_petsc_mat_base.tpp"
