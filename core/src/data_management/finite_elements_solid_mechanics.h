#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>
#include <tuple>

#include "data_management/data.h"
#include "control/types.h"
#include "mesh/mesh.h"
#include "field_variable/field_variable.h"
#include "equation/type_traits.h"

namespace Data
{

/** Common base class for solid mechanics finite elements, 1.) penalty formulation, 2.) mixed formulation
 */
template<typename BasisOnMeshType,typename Term>
class FiniteElementsSolidMechanics :
  public Data<BasisOnMeshType>
{
public:

  //! constructor
  FiniteElementsSolidMechanics(DihuContext context);

  //! destructor
  ~FiniteElementsSolidMechanics();

  //! initialize the object, create all stored data
  virtual void initialize() override;

  //! return reference to a stiffness matrix
  std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> tangentStiffnessMatrix();

  //! return reference to the residual field, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()> &residual();

  //! return reference to the increment field, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()> &increment(){LOG(FATAL) << "this should not be in use";}
  FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()> &solution(){LOG(FATAL) << "this should not be in use";}

  //! return reference to the actual geometry field, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<BasisOnMeshType,3> &geometryActual();

  //! return reference to the reference configuration geometry field, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<BasisOnMeshType,3> &geometryReference();

  //! return reference to the displacements field, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()> &displacements();

  //! return reference to the wExt field, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()> &externalVirtualWork();

  //! return reference to the wInt field, the PETSc Vec can be obtained via fieldVariable.values()
  FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()> &internalVirtualWork();

  //! alias for externalVirtualWork, needed such that rhs setting functionality works
  FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()> &rightHandSide();

  //! 3D vector to be used in 2D problems for adding to actual geometry which is also 3D
  Vec &fullIncrement();

  //! get a reference to the tangent stiffness matrix as it is used for the nonlinear solver
  std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> solverMatrixTangentStiffness();

  //! get a reference to the tangent stiffness matrix that is used for the finite difference approximation of the stiffness matrix if analytic jacobian is used
  std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> solverMatrixTangentStiffnessFiniteDifferences();

  //! get a reference to the residual vector as it is used for the nonlinear solver
  Vec &solverVariableResidual();

  //! get a reference to the solution vector as it is used for the nonlinear solver
  Vec &solverVariableSolution();

  //! get a reference to the reduced virtual work vector that contains entries of virtual work for indices that have no Dirichlet BC
  Vec &internalVirtualWorkReduced();

  //! perform the final assembly of petsc
  void finalAssembly();

  //! print all stored data to stdout
  void print();

  //! if the discretization matrix is already initialized
  bool massMatrixInitialized();

  //! create PETSc matrix
  void initializeMassMatrix();

  //! if the external virtual energy does not depend on deformation and thus is constant
  bool &externalVirtualWorkIsConstant();

  //! return a reference to the discretization matrix
  std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> massMatrix();

  //! get the number of unknows in the solution variable which is 3*nNodes
  virtual dof_no_t nUnknownsLocalWithGhosts();

  //! field variables that will be output by outputWriters
  typedef std::tuple<
    std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>>,  // geometryReference
    std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>>,  // actual geometry (stored in mesh)
    std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()>>,  // displacements
    std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()>>,   // residual
    std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()>>   // externalVirtualWork
  > OutputFieldVariables;

  //! get pointers to all field variables that can be written by output writers
  OutputFieldVariables getOutputFieldVariables();

  //! return reference to a stiffness matrix. This method is usually called for solving the linear system, but in this case we have an nonlinear system that does not retrieve the stiffness matrix
  std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> stiffnessMatrix(){LOG(FATAL) << "this should not be in use";}

  //! return the value of computeWithReducedVectors. If the vector of unknowns only contains the real degrees of freedom and not the variables with Dirichlet BCs. This is maybe slower because copying of data is required, but the system to solve is smaller
  bool computeWithReducedVectors();

  //! initialize the solver variables, solverMatrixTangentStiffness_, solverVariableSolution_ and solverVariableResidual_. @param nDofs the number of unknowns in the solver, i.e. reduced by Dirichlet BC
  void initializeSolverVariables(dof_no_t nDofs);

protected:

  //! initializes the vector and stiffness matrix with size
  void createPetscObjects();

  //! get maximum number of expected non-zeros in stiffness matrix
  void getPetscMemoryParameters(int &diagonalNonZeros, int &offdiagonalNonZeros);

  //! get the number of rows and columns to be used for setup of tangent stiffness matrix. This is different for mixed formulation.
  virtual const dof_no_t getTangentStiffnessMatrixNRows();

  PartitionedPetscMat<BasisOnMeshType> tangentStiffnessMatrix_;     ///< the tangent stiffness matrix which is the jacobian for the nonlinear problem. This is the non-reduced matrix that also contains entries for Dirichlet BC values. If computeWithReducedVectors_ is has to be reduced before computation.
  PartitionedPetscMat<BasisOnMeshType> solverMatrixTangentStiffnessFiniteDifferences_;    ///< this tangent stiffness matrix is used for the finite difference jacobian when an analytic jacobian gets computed in tangentStiffnessMatrix_
  PartitionedPetscMat<BasisOnMeshType> massMatrix_;  ///< a matrix that, applied to a rhs vector f, gives the rhs vector in weak formulation

  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()>> residual_;           ///< the residual vector, needed in the solution process by PETSc
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()>> increment_;           ///< the increments vector

  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>> geometryReference_;   //< geometry field in reference configuration, the geometry in actual configuration is stored by mesh_
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()>> displacements_;        //< current displacements
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()>> externalVirtualWork_;        //< the external virtual work vector δW_ext
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()>> internalVirtualWork_;        //< the internal virtual work vector δW_int

  Vec fullIncrement_;   ///< only for 2D problems this vec is a 3D vector that is filled from the 2D displacements vector and afterwards added to the geometry values.

  PartitionedPetscMat<BasisOnMeshType> solverMatrixTangentStiffness_;    ///< the tangent stiffness matrix used as jacobian in the nonlinear solver, in reduced form if we use reduce quantities (displacements without Dirichlet BC)
  Vec solverVariableResidual_;  ///< this vector is used to store a reduced version of the residual for the nonlinear solver
  Vec solverVariableSolution_;  ///< this vector is used to store the reduced displacements vector that does contain the same as displacements_ but without values for Dirichlet BCs
  Vec internalVirtualWorkReduced_;  ///< reduced vector for the internal virtual work

  bool externalVirtualWorkIsConstant_;   //< if the external virtual energy does not depend on deformation and thus is constant

  bool massMatrixInitialized_ = false;    ///< if the discretization matrix was initialized
  const bool computeWithReducedVectors_;    ///< if the vector of unknowns only contains the real degrees of freedom and not the variables with Dirichlet BCs. This is maybe slower because copying of data is required, but the system to solve is smaller.

};

/** inherit everything from FiniteElementsSolidMechanics
 */
template<typename BasisOnMeshType,typename Term>
class FiniteElements<
  BasisOnMeshType,
  Term,
  Equation::isSolidMechanics<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
> :
  public FiniteElementsSolidMechanics<BasisOnMeshType,Term>
{
public:
  //! inherit constructor
  using FiniteElementsSolidMechanics<BasisOnMeshType,Term>::FiniteElementsSolidMechanics;
};

}  // namespace

#include "data_management/finite_elements_solid_mechanics.tpp"
