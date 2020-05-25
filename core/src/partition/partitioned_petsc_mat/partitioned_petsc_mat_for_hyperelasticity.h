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
 *
 *  The template parameter nDisplacementComponents specifies the number of displacement components and therefore equals 3.
 *
 *  For incompressible formulation, where p is an unknown:
 *    The matrix has the following structure if nDisplacementComponents == 3:
 *
 *    (UU UP)   | 3 rows per dof
 *    (PU PP)   | 1 row per dof
 *
 *    The matrix has the following structure if nDisplacementComponents == 6:
 *
 *    (UU UV UP)   | 3 rows per dof
 *    (VU VV VP)   | 3 rows per dof
 *    (PU PV PP)   | 1 row per dof
 *
 * For compressible formulation:
 *    The matrix has the following structure if nDisplacementComponents == 3:
 *
 *    (UU)   | 3 rows per dof
 *
 *    The matrix has the following structure if nDisplacementComponents == 6:
 *
 *    (UU UV)   | 3 rows per dof
 *    (VU VV)   | 3 rows per dof
 *
 */

template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term, int nDisplacementComponents = 3>
class PartitionedPetscMatForHyperelasticityBase :
  public PartitionedPetscMatOneComponent<FunctionSpace::Generic, FunctionSpace::Generic>
{
public:

  //! constructor, create square sparse matrix, the number of entries is given by partitionedPetscVecForHyperelasticity
  PartitionedPetscMatForHyperelasticityBase(
    std::shared_ptr<PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nDisplacementComponents>> partitionedPetscVecForHyperelasticity,
    int nNonZerosDiagonal, int nNonZerosOffdiagonal,
    std::string name);

  //! this is the only special set function to set entries in the jacobian matrix (apart from the vectorized version, below).
  void setValue(int componentNoRow, PetscInt row, int componentNoColumn, PetscInt column, PetscScalar value, InsertMode mode);

  //! this is the only special set function to set entries in the jacobian matrix (apart from the non-vectorized version, above).
  //! Set the given values to all rows and columns of the respective components
  void setValue(int componentNoRow, Vc::int_v row, int componentNoColumn, Vc::int_v column, PetscScalar value, InsertMode mode);

  //! wrapper of MatSetValues for a vectorized value, sets matrix[rows[i],columns[i]] = values[i] for i = 1,..,Vc::double_v::size(), i.e. not matrix[rows[i],columns[j]] = ...
  void setValue(int componentNoRow, Vc::int_v rows, int componentNoColumn, Vc::int_v columns, Vc::double_v values, InsertMode mode);

  //! output the matrix to the file in globalNatural ordering, it has the same form regardless of number of ranks and therefore can be used to compare output with different ranks
  void dumpMatrixGlobalNatural(std::string filename);

  //! get a submatrix of this matrix. If nDisplacementComponents == 6, rowVariableNo and columnVariableNo can be 0,1 or 2 to select U,V or P.
  //! If the matrix does not contain velocities and nDisplacementComponents == 3, rowVariableNo and columnVariableNo can be 0 or 1.
  //! Then, getSubmatrix(0,0) = UU, getSubmatrix(0,1) = UP, getSubmatrix(1,0) = PU, getSubmatrix(1,1) = PP
  Mat getSubmatrix(int rowVariableNo, int columnVariableNo);

protected:

  std::shared_ptr<PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nDisplacementComponents>>
    partitionedPetscVecForHyperelasticity_;     //< one instance of the corresponding PartitionedPetscVecForHyperelasticity class, which handles all numbering
};

/** Class that adds possibility to dump matrix
 *  This is only possible for structured meshes, not for composite structured meshes
 */
template<typename DisplacementsFunctionSpaceType, typename PressureFunctionSpaceType, typename Term, int nDisplacementComponents = 3>
class PartitionedPetscMatForHyperelasticity :
  public PartitionedPetscMatForHyperelasticityBase<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nDisplacementComponents>
{
public:
  using PartitionedPetscMatForHyperelasticityBase<DisplacementsFunctionSpaceType,PressureFunctionSpaceType,Term,nDisplacementComponents>::PartitionedPetscMatForHyperelasticityBase;

  //! output the matrix to the file in globalNatural ordering, it has the same form regardless of number of ranks and therefore can be used to compare output with different ranks
  void dumpMatrixGlobalNatural(std::string filename){}
};

/** class that adds possibility to dump matrix, only for structured meshes, not composite meshes
 */
template<typename PressureFunctionSpaceType, typename Term, int nDisplacementComponents>
class PartitionedPetscMatForHyperelasticity<
  FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<2>>,
  PressureFunctionSpaceType,
  Term,
  nDisplacementComponents
> :
  public PartitionedPetscMatForHyperelasticityBase<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<2>>, PressureFunctionSpaceType, Term, nDisplacementComponents>
{
public:
  typedef FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<2>> DisplacementsFunctionSpaceType;

  using PartitionedPetscMatForHyperelasticityBase<FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>,BasisFunction::LagrangeOfOrder<2>>, PressureFunctionSpaceType, Term, nDisplacementComponents>::PartitionedPetscMatForHyperelasticityBase;

  //! output the matrix to the file in globalNatural ordering, it has the same form regardless of number of ranks and therefore can be used to compare output with different ranks
  void dumpMatrixGlobalNatural(std::string filename);
};

#include "partition/partitioned_petsc_mat/partitioned_petsc_mat_for_hyperelasticity.tpp"
