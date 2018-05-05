#include "data_management/finite_elements.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <memory>

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "control/dihu_context.h"
#include "utility/petsc_utility.h"
#include "basis_on_mesh/basis_on_mesh.h"
#include "mesh/unstructured_deformable.h"
#include "basis_function/hermite.h"

namespace Data
{

template<typename BasisOnMeshType,typename Term>
FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
FiniteElementsSolidMechanics(DihuContext context) : 
  Data<BasisOnMeshType>(context), computeWithReducedVectors_(true)  // true: works for analytic and numeric Jacobian in nonlinear solver, false: works only for analyitc Jacobian
{
  LOG(TRACE) << "Data::FiniteElements constructor";
  //PythonUtility::printDict(this->context_.getPythonConfig());
}

template<typename BasisOnMeshType,typename Term>
FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
~FiniteElementsSolidMechanics()
{
  PetscErrorCode ierr;
  // free PETSc objects
  if (this->initialized_)
  {
    ierr = MatDestroy(&this->tangentStiffnessMatrix_); CHKERRV(ierr);
  }
}

template<typename BasisOnMeshType,typename Term>
void FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
initialize()
{
  Data<BasisOnMeshType>::initialize();
  externalVirtualWorkIsConstant_ = true;
}

template<typename BasisOnMeshType,typename Term>
void FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
getPetscMemoryParameters(int &diagonalNonZeros, int &offdiagonalNonZeros)
{
  const int D = this->mesh_->dimension();
  const int nDofsPerNode = BasisOnMesh::BasisOnMeshBaseDim<1,typename BasisOnMeshType::BasisFunction>::nDofsPerNode();
  const int nDofsPerBasis = BasisOnMesh::BasisOnMeshBaseDim<1,typename BasisOnMeshType::BasisFunction>::nDofsPerElement();
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

template<typename BasisOnMeshType,typename Term>
dof_no_t FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
nUnknowns()
{
  return this->mesh_->nNodes() * BasisOnMeshType::dim();  // D components for displacements
}

template<typename BasisOnMeshType,typename Term>
const dof_no_t FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
getTangentStiffnessMatrixNRows()
{
  const int D = BasisOnMeshType::dim();
  return this->mesh_->nDofs() * D;
}

template<typename BasisOnMeshType,typename Term>
void FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
createPetscObjects()
{
  const int D = BasisOnMeshType::dim();
 
  this->residual_ = this->mesh_->template createFieldVariable<D>("residual");
  this->externalVirtualWork_ = this->mesh_->template createFieldVariable<D>("externalVirtualWork");
  this->internalVirtualWork_ = this->mesh_->template createFieldVariable<D>("internalVirtualWork");
  this->increment_ = this->mesh_->template createFieldVariable<D>("increment");
  this->displacements_ = this->mesh_->template createFieldVariable<D>("displacements");
  this->geometryReference_ = this->mesh_->template createFieldVariable<3>("geometryReference");
  
  // set geometryReference to be the same as the initial geometry field
  this->geometryReference_->setValues(this->mesh_->geometryField());
  
  // create a 3D vector
  PetscErrorCode ierr;
  if (BasisOnMeshType::dim() == 2)
  {
    ierr = VecDuplicate(this->mesh_->geometryField().values(), &fullIncrement_); CHKERRV(ierr);
  }
  
  // create PETSc matrix object
  
  // dimension of the tangent stiffness matrix
  dof_no_t tangentStiffnessMatrixNRows = this->getTangentStiffnessMatrixNRows();
  
  LOG(DEBUG)<<"FiniteElementsSolidMechanics<BasisOnMeshType,Term>::createPetscObjects, "
    << "dimension of tangent stiffness matrix: "<<tangentStiffnessMatrixNRows<<"x"<<tangentStiffnessMatrixNRows<<"";
  
  // PETSc MatCreateAIJ parameters
  int diagonalNonZeros = 3;   // number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
  int offdiagonalNonZeros = 0;   //  number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows)
  
  getPetscMemoryParameters(diagonalNonZeros, offdiagonalNonZeros);
  diagonalNonZeros = std::min(diagonalNonZeros, tangentStiffnessMatrixNRows);
  offdiagonalNonZeros = std::min(offdiagonalNonZeros, tangentStiffnessMatrixNRows);
  
  LOG(DEBUG) << "d="<<this->mesh_->dimension()
    <<", number of diagonal non-zeros: "<<diagonalNonZeros<<", number of off-diagonal non-zeros: "<<offdiagonalNonZeros;
  
  ierr = MatCreate(PETSC_COMM_WORLD, &this->tangentStiffnessMatrix_); CHKERRV(ierr);
  ierr = MatSetSizes(this->tangentStiffnessMatrix_, PETSC_DECIDE, PETSC_DECIDE, tangentStiffnessMatrixNRows, tangentStiffnessMatrixNRows); CHKERRV(ierr);
  ierr = MatSetFromOptions(this->tangentStiffnessMatrix_); CHKERRV(ierr);


  // dense matrix
  //ierr = MatSetUp(this->tangentStiffnessMatrix_); CHKERRV(ierr);
  
  // sparse matrix
  ierr = MatMPIAIJSetPreallocation(this->tangentStiffnessMatrix_, diagonalNonZeros, NULL, offdiagonalNonZeros, NULL); CHKERRV(ierr);
  ierr = MatSeqAIJSetPreallocation(this->tangentStiffnessMatrix_, diagonalNonZeros, NULL); CHKERRV(ierr);
  
  // allow additional non-zero entries in the stiffness matrix for UnstructuredDeformable mesh
  ierr = MatSetOption(this->tangentStiffnessMatrix_, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE); CHKERRV(ierr);
}

template<typename BasisOnMeshType,typename Term>
void FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
finalAssembly()
{
  PetscErrorCode ierr;
  // communicate portions to the right processors before using the matrix and vector in computations
  ierr = MatAssemblyBegin(this->tangentStiffnessMatrix_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  
  ierr = MatAssemblyEnd(this->tangentStiffnessMatrix_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  
  /*
  if (computeWithReducedVectors_)
  {
    ierr = MatAssemblyBegin(this->tangentStiffnessMatrixReduced_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
    ierr = MatAssemblyEnd(this->tangentStiffnessMatrixReduced_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  }
  */
  
  LOG(DEBUG) << "finalAssembly";
}

template<typename BasisOnMeshType,typename Term>
Mat &FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
tangentStiffnessMatrix()
{
  return this->tangentStiffnessMatrix_;
}

template<typename BasisOnMeshType,typename Term>
Mat &FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
solverMatrixTangentStiffnessFiniteDifferences()
{
  return this->solverMatrixTangentStiffnessFiniteDifferences_;
}

template<typename BasisOnMeshType,typename Term>
FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()> &FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
residual()
{
  return *this->residual_;
}

template<typename BasisOnMeshType,typename Term>
FieldVariable::FieldVariable<BasisOnMeshType,3> &FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
geometryActual()
{
  return this->mesh_->geometryField();
}

template<typename BasisOnMeshType,typename Term>
FieldVariable::FieldVariable<BasisOnMeshType,3> &FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
geometryReference()
{
  return *this->geometryReference_;
}

template<typename BasisOnMeshType,typename Term>
FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()> &FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
displacements()
{
  return *this->displacements_;
}

template<typename BasisOnMeshType,typename Term>
FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()> &FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
rightHandSide()
{
  return *this->externalVirtualWork_;
}

template<typename BasisOnMeshType,typename Term>
FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()> &FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
externalVirtualWork()
{
  return *this->externalVirtualWork_;
}

template<typename BasisOnMeshType,typename Term>
FieldVariable::FieldVariable<BasisOnMeshType,BasisOnMeshType::dim()> &FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
internalVirtualWork()
{
  return *this->internalVirtualWork_;
}

template<typename BasisOnMeshType,typename Term>
Mat &FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
massMatrix()
{
  return this->massMatrix_;
}

template<typename BasisOnMeshType,typename Term>
void FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
print()
{
  if (!VLOG_IS_ON(4))
    return;
  
  VLOG(4)<<"======================";
  int nRows, nColumns;
  MatGetSize(this->tangentStiffnessMatrix_, &nRows, &nColumns);
  VLOG(4)<<"tangentStiffnessMatrix ("<<nRows<<" x "<<nColumns<<") and rhs:";
  
  VLOG(4) << std::endl<<PetscUtility::getStringMatrixVector(this->tangentStiffnessMatrix_, this->residual_->values());
  VLOG(4) << "sparsity pattern: " << std::endl << PetscUtility::getStringSparsityPattern(this->tangentStiffnessMatrix_);
  
  MatInfo info;
  MatGetInfo(this->tangentStiffnessMatrix_, MAT_LOCAL, &info);
  
  VLOG(4)<<"Matrix info: "<<std::endl
    <<"block_size: "<<info.block_size<<std::endl
    <<"number of nonzeros: allocated: "<<info.nz_allocated<<", used: "<<info.nz_used<<", unneeded: "<<info.nz_unneeded<<std::endl
    <<"memory allocated: "<<info.memory<<std::endl
    <<"number of matrix assemblies called: "<<info.assemblies<<std::endl
    <<"number of mallocs during MatSetValues(): "<<info.mallocs<<std::endl
    <<"fill ratio for LU/ILU: given: "<<info.fill_ratio_given<<", needed: "<<info.fill_ratio_needed<<std::endl 
    <<"number of mallocs during factorization: "<<info.factor_mallocs<<std::endl;
    
    
  VLOG(4)<<"======================";
  
  int nEntries;
  VecGetSize(this->residual_->values(), &nEntries);
  VLOG(4)<<"residual ("<<nEntries<<" entries):";
  VLOG(4)<<PetscUtility::getStringVector(this->residual_->values());
  VLOG(4)<<"======================";
  
  VecGetSize(this->displacements_->values(), &nEntries);
  VLOG(4)<<"displacements ("<<nEntries<<" entries):";
  VLOG(4)<<PetscUtility::getStringVector(this->displacements_->values());
  VLOG(4)<<"======================";
}

template<typename BasisOnMeshType,typename Term>
bool FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
massMatrixInitialized()
{
  return this->massMatrixInitialized_;
}

template<typename BasisOnMeshType,typename Term>
bool &FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
externalVirtualWorkIsConstant()
{
  return this->externalVirtualWorkIsConstant_;
}

template<typename BasisOnMeshType,typename Term>
void FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
initializeMassMatrix()
{
  // determine problem size
  int nEntries;
  VecGetSize(this->solution_->values(), &nEntries);
  
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

template<typename BasisOnMeshType,typename Term>
typename FiniteElementsSolidMechanics<BasisOnMeshType,Term>::OutputFieldVariables 
FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
getOutputFieldVariables()
{
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>> actualGeometryField 
    = std::make_shared<FieldVariable::FieldVariable<BasisOnMeshType,3>>(this->mesh_->geometryField());
  /*
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>> generalField;
  
  generalField = std::static_pointer_cast<FieldVariable::FieldVariable<BasisOnMeshType,3>>(this->mesh_->fieldVariable("general"));
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


template<typename BasisOnMeshType,typename Term>
bool FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
computeWithReducedVectors()
{
  return computeWithReducedVectors_;
}
  
template<typename BasisOnMeshType,typename Term>
Vec &FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
fullIncrement()
{
  return fullIncrement_;
}

template<typename BasisOnMeshType,typename Term>
Mat &FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
solverMatrixTangentStiffness()
{
  return solverMatrixTangentStiffness_;
}

template<typename BasisOnMeshType,typename Term>
Vec &FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
solverVariableSolution()
{
  return solverVariableSolution_;
}

template<typename BasisOnMeshType,typename Term>
Vec &FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
solverVariableResidual()
{
  return solverVariableResidual_;
}

template<typename BasisOnMeshType,typename Term>
Vec &FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
internalVirtualWorkReduced()
{
  return internalVirtualWorkReduced_;
}

template<typename BasisOnMeshType,typename Term>
void FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
initializeSolverVariables(int nDofs)
{
  const int D = BasisOnMeshType::dim();
  const int nDofsDisplacementsReduced = this->mesh_->nDofs() * D;
 
  LOG(DEBUG) << "initializeSolverVariables with nDofs=" << nDofs << ", nDofsDisplacementsReduced=" << nDofsDisplacementsReduced;
  
  PetscErrorCode ierr;
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
