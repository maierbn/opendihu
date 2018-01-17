#include "data_management/finite_elements.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "control/dihu_context.h"
#include <utility/petsc_utility.h>

namespace Data
{

template<typename BasisOnMeshType>
FiniteElements<BasisOnMeshType>::
FiniteElements(const DihuContext &context) : Data<BasisOnMeshType>(context)
{
  this->disablePrinting_ = true;
  if (PythonUtility::containsKey(this->context_.getPythonConfig(), "disablePrinting"))
  {
    this->disablePrinting_ = PythonUtility::getOptionBool(this->context_.getPythonConfig(), "disablePrinting", false);
  }
  
  this->disableMatrixPrinting_ = true;
  if (PythonUtility::containsKey(this->context_.getPythonConfig(), "disableMatrixPrinting"))
  {
    this->disableMatrixPrinting_ = PythonUtility::getOptionBool(this->context_.getPythonConfig(), "disableMatrixPrinting", false);
  }
}

template<typename BasisOnMeshType>
FiniteElements<BasisOnMeshType>::
~FiniteElements()
{
  PetscErrorCode ierr;
  // free PETSc objects
  if (this->initialized_)
  {
    ierr = MatDestroy(&this->stiffnessMatrix_); CHKERRV(ierr);
  }
}

template<typename BasisOnMeshType>
void FiniteElements<BasisOnMeshType>::
createPetscObjects()
{
  element_idx_t n = this->mesh_->nNodes();
  
  LOG(DEBUG)<<"FiniteElements<BasisOnMeshType>::createPetscObjects("<<n<<")";
  
  this->rhs_ = this->mesh_->createFieldVariable("rhs");
  this->solution_ = this->mesh_->createFieldVariable("solution");
  
  PetscErrorCode ierr;
  // create PETSc matrix object
  
  // PETSc MatCreateAIJ parameters
  int d_nz = 3;   // number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
  int o_nz = 0;   //  number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows)
  
  switch (this->mesh_->dimension())
  {
  case 1:
    d_nz = 3;
    o_nz = 3;
    break;
  case 2:
    d_nz = 9;
    o_nz = 9;
    break;
  case 3:
    d_nz = 27;
    o_nz = 27;
    break;
  };
  
  LOG(DEBUG) << "d="<<this->mesh_->dimension()
    <<", number of diagonal non-zeros: "<<d_nz<<", number of off-diagonal non-zeros: "<<o_nz;
  
  ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, n, n, 
                      d_nz, NULL, o_nz, NULL, &this->stiffnessMatrix_); CHKERRV(ierr);
  ierr = MatMPIAIJSetPreallocation(this->stiffnessMatrix_, d_nz, NULL, o_nz, NULL); CHKERRV(ierr);
  
  // dense matrix:
  if (false)
  {
    ierr = MatCreate(PETSC_COMM_WORLD, &this->stiffnessMatrix_);  CHKERRV(ierr);
    ierr = MatSetSizes(this->stiffnessMatrix_, PETSC_DECIDE,PETSC_DECIDE, n, n);  CHKERRV(ierr);
    ierr = MatSetFromOptions(this->stiffnessMatrix_); CHKERRV(ierr);
    ierr = MatSetUp(this->stiffnessMatrix_); CHKERRV(ierr);
  }
}

template<typename BasisOnMeshType>
void FiniteElements<BasisOnMeshType>::
finalAssembly()
{
  PetscErrorCode ierr;
  // communicate portions to the right processors before using the matrix and vector in computations
  ierr = VecAssemblyBegin(this->rhs_->values()); CHKERRV(ierr);
  ierr = MatAssemblyBegin(this->stiffnessMatrix_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  
  ierr = VecAssemblyEnd(this->rhs_->values()); CHKERRV(ierr);
  ierr = MatAssemblyEnd(this->stiffnessMatrix_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  
  LOG(DEBUG) << "finalAssembly";
}

template<typename BasisOnMeshType>
Mat &FiniteElements<BasisOnMeshType>::
stiffnessMatrix()
{
  return this->stiffnessMatrix_;
}

template<typename BasisOnMeshType>
FieldVariable::FieldVariable<BasisOnMeshType> &FiniteElements<BasisOnMeshType>::
rightHandSide()
{
  return *this->rhs_;
}

template<typename BasisOnMeshType>
FieldVariable::FieldVariable<BasisOnMeshType> &FiniteElements<BasisOnMeshType>::
solution()
{
  return *this->solution_;
}

template<typename BasisOnMeshType>
Mat &FiniteElements<BasisOnMeshType>::
discretizationMatrix()
{
  return this->discretizationMatrix_;
}

template<typename BasisOnMeshType>
void FiniteElements<BasisOnMeshType>::
print()
{
  if (this->disablePrinting_)
    return;
  
  VLOG(1)<<"======================";
  int nRows, nColumns;
  MatGetSize(this->stiffnessMatrix_, &nRows, &nColumns);
  VLOG(1)<<"stiffnessMatrix ("<<nRows<<" x "<<nColumns<<") and rhs:";
  
  if (!this->disableMatrixPrinting_)
  {
    VLOG(1) << std::endl<<PetscUtility::getStringMatrixVector(this->stiffnessMatrix_, this->rhs_->values());
    VLOG(1) << "sparsity pattern: " << std::endl << PetscUtility::getStringSparsityPattern(this->stiffnessMatrix_);
  }
  
  MatInfo info;
  MatGetInfo(this->stiffnessMatrix_, MAT_LOCAL, &info);
  
  VLOG(1)<<"Matrix info: "<<std::endl
    <<"block_size: "<<info.block_size<<std::endl
    <<"number of nonzeros: allocated: "<<info.nz_allocated<<", used: "<<info.nz_used<<", unneeded: "<<info.nz_unneeded<<std::endl
    <<"memory allocated: "<<info.memory<<std::endl
    <<"number of matrix assemblies called: "<<info.assemblies<<std::endl
    <<"number of mallocs during MatSetValues(): "<<info.mallocs<<std::endl
    <<"fill ratio for LU/ILU: given: "<<info.fill_ratio_given<<", needed: "<<info.fill_ratio_needed<<std::endl 
    <<"number of mallocs during factorization: "<<info.factor_mallocs<<std::endl;
    
    
  VLOG(1)<<"======================";
  
  int nEntries;
  VecGetSize(this->rhs_->values(), &nEntries);
  VLOG(1)<<"rhs ("<<nEntries<<" entries):";
  VLOG(1)<<PetscUtility::getStringVector(this->rhs_->values());
  VLOG(1)<<"======================";
  
  VecGetSize(this->solution_->values(), &nEntries);
  VLOG(1)<<"solution ("<<nEntries<<" entries):";
  VLOG(1)<<PetscUtility::getStringVector(this->solution_->values());
  VLOG(1)<<"======================";
}

template<typename BasisOnMeshType>
bool FiniteElements<BasisOnMeshType>::
discretizationMatrixInitialized()
{
  return this->discretizationMatrixInitialized_;
}

template<typename BasisOnMeshType>
void FiniteElements<BasisOnMeshType>::
initializeDiscretizationMatrix()
{
  // determine problem size
  int nEntries;
  VecGetSize(this->rhs_->values(), &nEntries);
  
  // create PETSc matrix object
  
  // PETSc MatCreateAIJ parameters
  int d_nz = 3;   // number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
  int o_nz = 0;   //  number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows)
  
  switch (this->mesh_->dimension())
  {
  case 1:
    d_nz = 3;
    o_nz = 3;
    break;
  case 2:
    d_nz = 9;
    o_nz = 9;
    break;
  case 3:
    d_nz = 27;
    o_nz = 27;
    break;
  };
  
  PetscErrorCode ierr;
  ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, nEntries, nEntries, 
                      d_nz, NULL, o_nz, NULL, &this->discretizationMatrix_); CHKERRV(ierr);
  ierr = MatMPIAIJSetPreallocation(this->discretizationMatrix_, d_nz, NULL, o_nz, NULL); CHKERRV(ierr);
  
  this->discretizationMatrixInitialized_ = true;
}

} // namespace Data
