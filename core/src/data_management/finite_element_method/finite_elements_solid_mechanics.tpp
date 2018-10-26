#include "data_management/finite_element_method/finite_elements.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <memory>

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "control/dihu_context.h"
#include "utility/petsc_utility.h"
#include "function_space/function_space.h"
#include "mesh/unstructured_deformable.h"
#include "basis_function/hermite.h"

namespace Data
{

template<typename FunctionSpaceType,typename Term>
FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
FiniteElementsSolidMechanics(DihuContext context) :
  Data<FunctionSpaceType>(context), computeWithReducedVectors_(true)  // true: works for analytic and numeric Jacobian in nonlinear solver, false: works only for analyitc Jacobian
{
  LOG(TRACE) << "Data::FiniteElements constructor";
  //PythonUtility::printDict(this->context_.getPythonConfig());
}

template<typename FunctionSpaceType,typename Term>
FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
~FiniteElementsSolidMechanics()
{
  PetscErrorCode ierr;
  // free PETSc objects
  if (this->initialized_)
  {
    ierr = MatDestroy(&this->tangentStiffnessMatrix_); CHKERRV(ierr);
  }
}

template<typename FunctionSpaceType,typename Term>
void FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
initialize()
{
  Data<FunctionSpaceType>::initialize();
  externalVirtualWorkIsConstant_ = true;
}

template<typename FunctionSpaceType,typename Term>
void FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
getPetscMemoryParameters(int &diagonalNonZeros, int &offdiagonalNonZeros)
{
  const int D = this->functionSpace_->dimension();
  const int nDofsPerNode = FunctionSpace::FunctionSpaceBaseDim<1,typename FunctionSpaceType::BasisFunction>::nDofsPerNode();
  const int nDofsPerBasis = FunctionSpace::FunctionSpaceBaseDim<1,typename FunctionSpaceType::BasisFunction>::nDofsPerElement();
  const int nOverlaps = (nDofsPerBasis*2 - 1) * nDofsPerNode;   // number of nodes of 2 neighbouring 1D elements (=number of ansatz functions in support of center ansatz function)

  // due to PETSc storage diagonalNonZeros and offdiagonalNonZeros should be both set to the maximum number of non-zero entries per row

  switch (D)
  {
  case 1:
    diagonalNonZeros = nOverlaps;
    offdiagonalNonZeros = diagonalNonZeros;
    break;
  case 2:
    diagonalNonZeros = pow(nOverlaps, 2) + 16;   // because of boundary conditions there can be more entries, which are all zero, but stored as non-zero
    offdiagonalNonZeros = diagonalNonZeros;
    break;
  case 3:
    diagonalNonZeros = pow(nOverlaps, 3);
    offdiagonalNonZeros = diagonalNonZeros;
    break;
  };
}

template<typename FunctionSpaceType,typename Term>
dof_no_t FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
nUnknownsLocalWithGhosts()
{
  return this->functionSpace_->nNodesLocalWithGhosts() * FunctionSpaceType::dim();  // D components for displacements
}

template<typename FunctionSpaceType,typename Term>
dof_no_t FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
getTangentStiffnessMatrixNRows()
{
  const int D = FunctionSpaceType::dim();
  return this->functionSpace_->nDofsLocal() * D;
}

template<typename FunctionSpaceType,typename Term>
void FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
createPetscObjects()
{
  const int D = FunctionSpaceType::dim();

  this->residual_ = this->functionSpace_->template createFieldVariable<D>("residual");
  this->externalVirtualWork_ = this->functionSpace_->template createFieldVariable<D>("externalVirtualWork");
  this->internalVirtualWork_ = this->functionSpace_->template createFieldVariable<D>("internalVirtualWork");
  this->increment_ = this->functionSpace_->template createFieldVariable<D>("increment");
  this->displacements_ = this->functionSpace_->template createFieldVariable<D>("displacements");
  this->geometryReference_ = this->functionSpace_->template createFieldVariable<3>("geometryReference");

  // set geometryReference to be the same as the initial geometry field
  this->geometryReference_->setValues(this->functionSpace_->geometryField());

  // create a 3D vector
  PetscErrorCode ierr;
  if (FunctionSpaceType::dim() == 2)
  {
    ierr = VecDuplicate(this->functionSpace_->geometryField().valuesGlobal(), &fullIncrement_); CHKERRV(ierr);
  }

  // create PETSc matrix object

  // dimension of the tangent stiffness matrix
  dof_no_t tangentStiffnessMatrixNRows = this->getTangentStiffnessMatrixNRows();

  LOG(DEBUG) << "FiniteElementsSolidMechanics<FunctionSpaceType,Term>::createPetscObjects, "
    << "dimension of tangent stiffness matrix: " << tangentStiffnessMatrixNRows << "x" << tangentStiffnessMatrixNRows << "";

  // PETSc MatCreateAIJ parameters
  int diagonalNonZeros = 3;   // number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
  int offdiagonalNonZeros = 0;   //  number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows)

  getPetscMemoryParameters(diagonalNonZeros, offdiagonalNonZeros);
  diagonalNonZeros = std::min(diagonalNonZeros, tangentStiffnessMatrixNRows);
  offdiagonalNonZeros = std::min(offdiagonalNonZeros, tangentStiffnessMatrixNRows);

  LOG(DEBUG) << "d=" << this->functionSpace_->dimension()
    << ", number of diagonal non-zeros: " << diagonalNonZeros << ", number of off-diagonal non-zeros: " <<offdiagonalNonZeros; 
  const int dimension = FunctionSpaceType::dim();
    
  const int nComponents = dimension;
  tangentStiffnessMatrix_ = std::make_shared<PartitionedPetscMat>(this->functionSpace_->meshPartition(), nComponents, diagonalNonZeros, offdiagonalNonZeros);
  
  // allow additional non-zero entries in the stiffness matrix for UnstructuredDeformable mesh
  ierr = MatSetOption(this->tangentStiffnessMatrix_.valuesGlobal(), MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE); CHKERRV(ierr);
}

template<typename FunctionSpaceType,typename Term>
void FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
finalAssembly()
{
  // communicate portions to the right processors before using the matrix and vector in computations
  this->tangentStiffnessMatrix_.assembly(MAT_FINAL_ASSEMBLY);
 
  /*
  if (computeWithReducedVectors_)
  {
    ierr = MatAssemblyBegin(this->tangentStiffnessMatrixReduced_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
    ierr = MatAssemblyEnd(this->tangentStiffnessMatrixReduced_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  }
  */

  LOG(DEBUG) << "finalAssembly";
}

template<typename FunctionSpaceType,typename Term>
std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
tangentStiffnessMatrix()
{
  return this->tangentStiffnessMatrix_;
}

template<typename FunctionSpaceType,typename Term>
std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
solverMatrixTangentStiffnessFiniteDifferences()
{
  return this->solverMatrixTangentStiffnessFiniteDifferences_;
}

template<typename FunctionSpaceType,typename Term>
FieldVariable::FieldVariable<FunctionSpaceType,FunctionSpaceType::dim()> &FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
residual()
{
  return *this->residual_;
}

template<typename FunctionSpaceType,typename Term>
FieldVariable::FieldVariable<FunctionSpaceType,3> &FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
geometryActual()
{
  return this->functionSpace_->geometryField();
}

template<typename FunctionSpaceType,typename Term>
FieldVariable::FieldVariable<FunctionSpaceType,3> &FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
geometryReference()
{
  return *this->geometryReference_;
}

template<typename FunctionSpaceType,typename Term>
FieldVariable::FieldVariable<FunctionSpaceType,FunctionSpaceType::dim()> &FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
displacements()
{
  return *this->displacements_;
}

template<typename FunctionSpaceType,typename Term>
FieldVariable::FieldVariable<FunctionSpaceType,FunctionSpaceType::dim()> &FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
rightHandSide()
{
  return *this->externalVirtualWork_;
}

template<typename FunctionSpaceType,typename Term>
FieldVariable::FieldVariable<FunctionSpaceType,FunctionSpaceType::dim()> &FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
externalVirtualWork()
{
  return *this->externalVirtualWork_;
}

template<typename FunctionSpaceType,typename Term>
FieldVariable::FieldVariable<FunctionSpaceType,FunctionSpaceType::dim()> &FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
internalVirtualWork()
{
  return *this->internalVirtualWork_;
}

template<typename FunctionSpaceType,typename Term>
std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
massMatrix()
{
  return this->massMatrix_;
}

template<typename FunctionSpaceType,typename Term>
void FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
print()
{
  if (!VLOG_IS_ON(4))
    return;

  VLOG(4) << "======================";
  int nRows, nColumns;
  MatGetSize(this->tangentStiffnessMatrix_, &nRows, &nColumns);
  VLOG(4) << "tangentStiffnessMatrix (" <<nRows << " x " <<nColumns << ") and rhs:";

  VLOG(4) << std::endl<<PetscUtility::getStringMatrixVector(this->tangentStiffnessMatrix_, this->residual_->values());
  VLOG(4) << "sparsity pattern: " << std::endl << PetscUtility::getStringSparsityPattern(this->tangentStiffnessMatrix_);

  MatInfo info;
  MatGetInfo(this->tangentStiffnessMatrix_, MAT_LOCAL, &info);

  VLOG(4) << "Matrix info: " << std::endl
    << "block_size: " <<info.block_size << std::endl
    << "number of nonzeros: allocated: " <<info.nz_allocated<< ", used: " <<info.nz_used<< ", unneeded: " <<info.nz_unneeded<< std::endl
    << "memory allocated: " <<info.memory<< std::endl
    << "number of matrix assemblies called: " <<info.assemblies << std::endl
    << "number of mallocs during MatSetValues(): " <<info.mallocs << std::endl
    << "fill ratio for LU/ILU: given: " <<info.fill_ratio_given<< ", needed: " <<info.fill_ratio_needed<< std::endl
    << "number of mallocs during factorization: " <<info.factor_mallocs << std::endl;


  VLOG(4) << "======================";

  int nEntries;
  VecGetSize(this->residual_->valuesLocal(), &nEntries);
  VLOG(4) << "residual (" <<nEntries << " entries):";
  VLOG(4) <<PetscUtility::getStringVector(this->residual_->valuesLocal());
  VLOG(4) << "======================";

  VecGetSize(this->displacements_->valuesLocal(), &nEntries);
  VLOG(4) << "displacements (" <<nEntries << " entries):";
  VLOG(4) <<PetscUtility::getStringVector(this->displacements_->valuesLocal());
  VLOG(4) << "======================";
}

template<typename FunctionSpaceType,typename Term>
bool FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
massMatrixInitialized()
{
  return this->massMatrixInitialized_;
}

template<typename FunctionSpaceType,typename Term>
bool &FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
externalVirtualWorkIsConstant()
{
  return this->externalVirtualWorkIsConstant_;
}

template<typename FunctionSpaceType,typename Term>
void FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
initializeMassMatrix()
{
  // determine problem size
  int nEntries;
  VecGetSize(this->solution_->valuesGlobal(), &nEntries);

  // create PETSc matrix object

  // PETSc MatCreateAIJ parameters
  int diagonalNonZeros = 3;   // number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
  int offdiagonalNonZeros = 0;   //  number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows)

  getPetscMemoryParameters(diagonalNonZeros, offdiagonalNonZeros);

  PetscErrorCode ierr;
  ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, nEntries, nEntries,
                      diagonalNonZeros, NULL, offdiagonalNonZeros, NULL, &this->massMatrix_); CHKERRV(ierr);
  ierr = MatMPIAIJSetPreallocation(this->massMatrix_, diagonalNonZeros, NULL, offdiagonalNonZeros, NULL); CHKERRV(ierr);

  this->massMatrixInitialized_ = true;
}

template<typename FunctionSpaceType,typename Term>
typename FiniteElementsSolidMechanics<FunctionSpaceType,Term>::OutputFieldVariables
FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
getOutputFieldVariables()
{
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> actualGeometryField
    = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,3>>(this->functionSpace_->geometryField());
  /*
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> generalField;

  generalField = std::static_pointer_cast<FieldVariable::FieldVariable<FunctionSpaceType,3>>(this->functionSpace_->fieldVariable("general"));
  if (!generalField)
   generalField = geometryField;
  */
  return OutputFieldVariables(
    geometryReference_,
    actualGeometryField,
    displacements_,
    residual_,
    externalVirtualWork_
  );
}


template<typename FunctionSpaceType,typename Term>
bool FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
computeWithReducedVectors()
{
  return computeWithReducedVectors_;
}

template<typename FunctionSpaceType,typename Term>
Vec &FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
fullIncrement()
{
  return fullIncrement_;
}

template<typename FunctionSpaceType,typename Term>
std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
solverMatrixTangentStiffness()
{
  return solverMatrixTangentStiffness_;
}

template<typename FunctionSpaceType,typename Term>
Vec &FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
solverVariableSolution()
{
  return solverVariableSolution_;
}

template<typename FunctionSpaceType,typename Term>
Vec &FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
solverVariableResidual()
{
  return solverVariableResidual_;
}

template<typename FunctionSpaceType,typename Term>
Vec &FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
internalVirtualWorkReduced()
{
  return internalVirtualWorkReduced_;
}

template<typename FunctionSpaceType,typename Term>
void FiniteElementsSolidMechanics<FunctionSpaceType,Term>::
initializeSolverVariables(int nDofs)
{
  const int D = FunctionSpaceType::dim();
  const int nDofsDisplacementsReduced = this->functionSpace_->nDofsLocal() * D;

  LOG(DEBUG) << "initializeSolverVariables with nDofs=" << nDofs << ", nDofsDisplacementsReduced=" << nDofsDisplacementsReduced;

  PetscErrorCode ierr;
  
  PetscUtility::createVector(this->solverVariableSolution_, nDofs, "solverVariableSolution", this->functionSpace_->partition());
  PetscUtility::createVector(this->solverVariableResidual_, nDofs, "solverVariableResidual", this->functionSpace_->partition());
  PetscUtility::createVector(this->internalVirtualWorkReduced_, nDofsDisplacementsReduced, "internalVirtualWorkReduced", this->functionSpace_->partition());
  /*
  ierr = VecCreate(PETSC_COMM_WORLD, &this->solverVariableSolution_);  CHKERRV(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD, &this->solverVariableResidual_);  CHKERRV(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD, &this->internalVirtualWorkReduced_);  CHKERRV(ierr);
  ierr = PetscObjectSetName((PetscObject) this->solverVariableSolution_, "solverVariableSolution"); CHKERRV(ierr);
  ierr = PetscObjectSetName((PetscObject) this->solverVariableResidual_, "solverVariableResidual"); CHKERRV(ierr);
  ierr = PetscObjectSetName((PetscObject) this->internalVirtualWorkReduced_, "internalVirtualWorkReduced"); CHKERRV(ierr);

  // initialize size of vectors
  ierr = VecSetSizes(this->solverVariableSolution_, PETSC_DECIDE, nDofs); CHKERRV(ierr);
  ierr = VecSetSizes(this->solverVariableResidual_, PETSC_DECIDE, nDofs); CHKERRV(ierr);
  ierr = VecSetSizes(this->internalVirtualWorkReduced_, PETSC_DECIDE, nDofsDisplacementsReduced); CHKERRV(ierr);

  // set sparsity type and other options
  ierr = VecSetFromOptions(this->solverVariableSolution_);  CHKERRV(ierr);
  ierr = VecSetFromOptions(this->solverVariableResidual_);  CHKERRV(ierr);
  ierr = VecSetFromOptions(this->internalVirtualWorkReduced_);  CHKERRV(ierr);
*/
  // set tangentStiffnessMatrixReduced_ to NULL, it is initialized the first time it is reduced from the full matrix using MatGetSubMatrix
  this->solverMatrixTangentStiffness_ = PETSC_NULL;

  /*
  // initialize reduced tangent stiffness matrix
  int diagonalNonZeros = 3;   // number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
  int offdiagonalNonZeros = 0;   //  number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows)

  getPetscMemoryParameters(diagonalNonZeros, offdiagonalNonZeros);
  diagonalNonZeros = std::min(diagonalNonZeros, nDofs);
  offdiagonalNonZeros = std::min(offdiagonalNonZeros, nDofs);

  ierr = MatCreate(PETSC_COMM_WORLD, &this->solverMatrixTangentStiffness_); CHKERRV(ierr);
  ierr = MatSetSizes(this->solverMatrixTangentStiffness_, PETSC_DECIDE, PETSC_DECIDE, nDofs, nDofs); CHKERRV(ierr);
  ierr = MatSetFromOptions(this->solverMatrixTangentStiffness_); CHKERRV(ierr);

  ierr = MatMPIAIJSetPreallocation(this->solverMatrixTangentStiffness_, diagonalNonZeros, NULL, offdiagonalNonZeros, NULL); CHKERRV(ierr);
  ierr = MatSeqAIJSetPreallocation(this->solverMatrixTangentStiffness_, diagonalNonZeros, NULL); CHKERRV(ierr);

  // allow additional non-zero entries in the stiffness matrix for UnstructuredDeformable mesh
  ierr = MatSetOption(this->solverMatrixTangentStiffness_, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE); CHKERRV(ierr);
  */
}

} // namespace Data
