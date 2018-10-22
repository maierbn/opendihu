#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>

#include "control/types.h"
#include "mesh/mesh.h"
#include "data_management/data.h"
#include "data_management/time_stepping.h"
#include "field_variable/field_variable.h"
#include "partition/partitioned_petsc_vec/partitioned_petsc_vec.h"
#include "partition/partitioned_petsc_mat/partitioned_petsc_mat.h"

namespace Data
{

/**  The datastructures used for implicit timestepping schemes.
  */
template<typename FunctionSpaceType, int nComponents>
class TimeSteppingImplicit : public TimeStepping<FunctionSpaceType,nComponents>
{
public:

  typedef FieldVariable::FieldVariable<FunctionSpaceType,nComponents> FieldVariableType;

  //! constructor of base class
  using TimeStepping<FunctionSpaceType,nComponents>::TimeStepping;

  //! return a reference to the rhs summand vector which is needed to apply the boundary conditions, the PETSc Vec can be obtained via fieldVariable->valuesGlobal()
  std::shared_ptr<FieldVariableType> boundaryConditionsRightHandSideSummand();

  //! return a reference to the rhs vector of the system, the PETSc Vec can be obtained via fieldVariable->valuesGlobal()
  std::shared_ptr<FieldVariableType> systemRightHandSide();

  //! get the system matrix required by the implicit time stepping
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> systemMatrix();

  //! get the integration matrix if required for the implicit time stepping
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> integrationMatrixRightHandSide();

  //! initialize the sytem matrix from a PETSc matrix that was already created, in this case by a MatMatMult
  void initializeSystemMatrix(Mat &systemMatrix);

  //! initialize the integration matrix on the rhs from a PETSc matrix that was already created, in this case by a MatConvert
  void initializeIntegrationMatrixRightHandSide(Mat &integrationMatrix);

  //! initializes a PETSc matrix that is already created, by other PETSc routines like MatConvert or MatMatMult
  void initializeMatrix(Mat &matrixIn, std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> matrixOut, std::string name);

  //! print all stored data to stdout
  void print() override;

private:

  //! initializes the vectors with size
  void createPetscObjects() override;

  std::shared_ptr<FieldVariableType> boundaryConditionsRightHandSideSummand_; ///< what needs to be added to rhs when dirichlet boundary conditions are applied
  std::shared_ptr<FieldVariableType> systemRightHandSide_;                    ///< the rhs vector of the implicit system

  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> systemMatrix_;  ///< the system matrix for implicit time stepping
  std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> integrationMatrixRightHandSide_; ///< the integration matrix for the right hand side if required by the implicit time stepping


};

} // namespace Data

#include "data_management/time_stepping_implicit.tpp"
