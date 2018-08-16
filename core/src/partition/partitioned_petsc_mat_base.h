#pragma once

#include <memory>

#include "control/types.h"
#include "partition/rank_subset.h"

/** Base class for a partitioned PetscMat
 */
class PartitionedPetscMatBase
{
public:
 
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
  
  //! get a reference to the PETSc matrix
  Mat &values();
  
protected:
 
  Mat matrix_;   ///< the global Petsc matrix, access it performed using MatSetValuesLocal() with local indices
};
