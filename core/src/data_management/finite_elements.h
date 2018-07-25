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
  ~FiniteElements();

  //! initialize the object, create all stored data
  virtual void initialize() override;

  //! return reference to a stiffness matrix
  std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> stiffnessMatrix();

  //! return reference to a right hand side vector, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<BasisOnMeshType,1> &rightHandSide();

  //! return reference to solution of the system, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<BasisOnMeshType,1> &solution();

  //! perform the final assembly of petsc
  void finalAssembly();

  //! print all stored data to stdout
  void print();

  //! if the discretization matrix is already initialized
  bool massMatrixInitialized();

  //! create PETSc matrix
  void initializeMassMatrix();

  //! return a reference to the discretization matrix
  std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> massMatrix();

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

  std::shared_ptr<PartitionedPetscMat> stiffnessMatrix_;     ///< the standard stiffness matrix of the finite element formulation
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,1>> rhs_;                 ///< the rhs vector in weak formulation
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,1>> solution_;            ///< the vector of the quantity of interest, e.g. displacement
  PartitionedPetscMat<BasisOnMeshType> massMatrix_;  ///< a matrix that, applied to a rhs vector f, gives the rhs vector in weak formulation

  bool massMatrixInitialized_ = false;    ///< if the discretization matrix was initialized

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
