#pragma once

#include <Python.h>  // has to be the first included header
#include <memory>

#include "control/types.h"
#include "partition/partitioned_petsc_mat/partitioned_petsc_mat_one_component.h"
#include "function_space/function_space_generic.h"
#include "function_space/function_space.h"
#include "partition/partitioned_petsc_vec/02_partitioned_petsc_vec_for_hyperelasticity.h"

/** A Petsc Mat to be used as stiffness matrix / jacobian in the nonlinear problem
 *  for the quasi_static_hyperelasticity_solver.
 *
 *  Everything related to indexing and number of elements is given by a vector of type
 *  PartitionedPetscVecForHyperelasticity. Read that class to understand how it is done.
 */

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType>
class PartitionedPetscMatForHyperelasticity :
  public PartitionedPetscMatOneComponent<FunctionSpace::Generic, FunctionSpace::Generic, typename FunctionSpace::Generic::Mesh>
{
public:

  //! constructor, create square sparse matrix, the number of entries is given by partitionedPetscVecForHyperelasticity
  PartitionedPetscMatForHyperelasticity(
    std::shared_ptr<PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType>> partitionedPetscVecForHyperelasticity,
    int diagonalNonZeros, int offdiagonalNonZeros,
    std::string name);

  //! this is the only special set function to set entries in the jacobian matrix.
  void setValue(int componentNoRow, PetscInt row, int componentNoColumn, PetscInt column, PetscScalar value, InsertMode mode);

  //! output the matrix to the file in globalNatural ordering, it has the same form regardless of number of ranks and therefore can be used to compare output with different ranks
  void dumpMatrixGlobalNatural(std::string filename);

  //! get a submatrix of the upper left part (only displacements)
  Mat getSubmatrixUU();

  //! get a submatrix of the upper right part
  Mat getSubmatrixUP();

  //! get a submatrix of the lower left part
  Mat getSubmatrixPU();

protected:

  std::shared_ptr<PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType>> partitionedPetscVecForHyperelasticity_;
};

#include "partition/partitioned_petsc_mat/partitioned_petsc_mat_for_hyperelasticity.tpp"
