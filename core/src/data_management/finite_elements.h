#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>
#include <tuple>

#include "data_management/data.h"
#include "data_management/diffusion_tensor.h"
#include "control/types.h"
#include "mesh/mesh.h"
#include "field_variable/field_variable.h"
#include "basis_on_mesh/mixed_basis_on_mesh.h"
#include "partition/partitioned_petsc_vec.h"
#include "partition/partitioned_petsc_mat.h"

namespace Data
{

template<typename BasisOnMeshType,typename Term,typename = Term,typename = typename BasisOnMeshType::BasisFunction>
class FiniteElements :
  public Data<BasisOnMeshType>,
  public DiffusionTensor<BasisOnMeshType::dim()>
{
public:

  //! constructor
  FiniteElements(DihuContext context);

  //! destructor
  virtual ~FiniteElements();

  //! initialize the object, create all stored data
  virtual void initialize() override;

  //! return reference to a right hand side vector, the PETSc Vec can be obtained via fieldVariable.valuesGlobal()
  FieldVariable::FieldVariable<BasisOnMeshType,1> &rightHandSide();

  //! return reference to solution of the system, the PETSc Vec can be obtained via fieldVariable.valuesGlobal()
  FieldVariable::FieldVariable<BasisOnMeshType,1> &solution();

  //! perform the final assembly of petsc
  void finalAssembly();

  //! print all stored data to stdout
  void print();
  
  //! create PETSc matrix
  void initializeMassMatrix();

  //! create the inverse of the lumped mass matrix
  void initializeInverseLumpedMassMatrix();

  //! initialize the sytem matrix from a PETSc matrix that was already created, in this case by a MatMatMult
  void initializeSystemMatrix(Mat &systemMatrix);

  //! return reference to a stiffness matrix
  std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> stiffnessMatrix();

  //! get the mass matrix
  std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> massMatrix();

  //! get the system matrix, corresponding to the specific time integration. (I - d*tM^(-1)*K) for the implicit Euler scheme.
  std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> systemMatrix();

  //! get the inversed lumped mass matrix
  std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> inverseLumpedMassMatrix();
  
  //! field variables that will be output by outputWriters
  typedef std::tuple<
    std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>>,  // geometry
    std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,1>>,  // solution
    std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,1>>   // rhs
  > OutputFieldVariables;

  //! get pointers to all field variables that can be written by output writers
  OutputFieldVariables getOutputFieldVariables();

private:

  //! initializes the vectors and stiffness matrix with size
  void createPetscObjects();

  //! get maximum number of expected non-zeros in stiffness matrix
  void getPetscMemoryParameters(int &diagonalNonZeros, int &offdiagonalNonZeros);

  std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> stiffnessMatrix_;      ///< the standard stiffness matrix of the finite element formulation
  std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> massMatrix_;           ///< the standard mass matrix, which is a matrix that, applied to a rhs vector f, gives the rhs vector in weak formulation
  std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> systemMatrix_;         ///< the system matrix for implicit time stepping, (I - dt*M^-1*K)
  std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> inverseLumpedMassMatrix_;         ///< the inverse lumped mass matrix that has only entries on the diagonal, they are the reciprocal of the row sums of the mass matrix

  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,1>> rhs_;                 ///< the rhs vector in weak formulation
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,1>> solution_;            ///< the vector of the quantity of interest, e.g. displacement

};

/*
#include "equation/type_traits.h"

template<typename BasisOnMeshType,typename Term>
class FiniteElements<
  BasisOnMeshType,
  Term,
  Equation::isSolidMechanics<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
> :
  public Data<BasisOnMeshType>
{
public:
  void tangentStiffnessMatrix();
};
*/
}  // namespace

#include "data_management/finite_elements.tpp"
#include "data_management/finite_elements_mixed.h"
#include "data_management/finite_elements_solid_mechanics.h"
