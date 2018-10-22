#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>
#include <tuple>

#include "data_management/data.h"
#include "control/types.h"
#include "mesh/mesh.h"
#include "field_variable/field_variable.h"
#include "function_space/mixed_function_space.h"
#include "partition/partitioned_petsc_vec/partitioned_petsc_vec.h"
#include "partition/partitioned_petsc_mat/partitioned_petsc_mat.h"

namespace Data
{

template<typename FunctionSpaceType>
class FiniteElementsBase :
  public Data<FunctionSpaceType>
{
public:

  //! constructor
  FiniteElementsBase(DihuContext context);

  //! destructor
  virtual ~FiniteElementsBase();

  //! initialize the object, create all stored data
  virtual void initialize() override;

  //! return reference to a right hand side vector, the PETSc Vec can be obtained via fieldVariable.valuesGlobal()
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> rightHandSide();

  //! return reference to solution of the system, the PETSc Vec can be obtained via fieldVariable.valuesGlobal()
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> solution();

  //! print all stored data to stdout
  void print();
  
  //! create PETSc matrix
  void initializeMassMatrix();

  //! create the inverse of the lumped mass matrix
  void initializeInverseLumpedMassMatrix();

  //! return reference to a stiffness matrix
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> stiffnessMatrix();

  //! get the mass matrix
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> massMatrix();

  //! get the inversed lumped mass matrix
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> inverseLumpedMassMatrix();
  
  //! field variables that will be output by outputWriters
  typedef std::tuple<
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>>,  // geometry
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>,  // solution
    std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>   // rhs
  > OutputFieldVariables;

  //! get pointers to all field variables that can be written by output writers
  OutputFieldVariables getOutputFieldVariables();

private:

  //! initializes the vectors and stiffness matrix with size
  void createPetscObjects();

  //! get maximum number of expected non-zeros in stiffness matrix
  void getPetscMemoryParameters(int &diagonalNonZeros, int &offdiagonalNonZeros);

  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> stiffnessMatrix_;      ///< the standard stiffness matrix of the finite element formulation
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> massMatrix_;           ///< the standard mass matrix, which is a matrix that, applied to a rhs vector f, gives the rhs vector in weak formulation
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> inverseLumpedMassMatrix_;         ///< the inverse lumped mass matrix that has only entries on the diagonal, they are the reciprocal of the row sums of the mass matrix

  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> rhs_;                 ///< the rhs vector in weak formulation
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> solution_;            ///< the vector of the quantity of interest, e.g. displacement

};

/*
#include "equation/type_traits.h"

template<typename FunctionSpaceType,typename Term>
class FiniteElements<
  FunctionSpaceType,
  Term,
  Equation::isSolidMechanics<Term>,
  BasisFunction::isNotMixed<typename FunctionSpaceType::BasisFunction>
> :
  public Data<FunctionSpaceType>
{
public:
  void tangentStiffnessMatrix();
};
*/
}  // namespace

#include "data_management/finite_elements_base.tpp"
