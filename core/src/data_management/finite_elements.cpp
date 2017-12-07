#include "data_management/finite_elements.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "easylogging++.h"

#include "control/python_utility.h"
#include "control/dihu_context.h"
#include <control/petsc_utility.h>

namespace Data
{

FiniteElements::FiniteElements(const DihuContext &context) : Data(context)
{
  disablePrinting_ = PythonUtility::getOptionBool(context_.getPythonConfig(), "disablePrinting", false);
  disableMatrixPrinting_ = PythonUtility::getOptionBool(context_.getPythonConfig(), "disableMatrixPrinting", false);
}

FiniteElements::~FiniteElements()
{
  PetscErrorCode ierr;
  // free PETSc objects
  if (initialized_)
  {
    ierr = MatDestroy(&stiffnessMatrix_); CHKERRV(ierr);
    ierr = VecDestroy(&rhs_); CHKERRV(ierr);
    ierr = VecDestroy(&solution_); CHKERRV(ierr);
  }
}

void FiniteElements::createPetscObjects()
{
  element_idx_t n = mesh_->nNodes();
  
  LOG(DEBUG)<<"FiniteElements::initVectorSize("<<n<<")"<<std::endl;
  
  PetscErrorCode ierr;
  // create PETSc vector object
  ierr = VecCreate(PETSC_COMM_WORLD, &rhs_);  CHKERRV(ierr);
  ierr = PetscObjectSetName((PetscObject) rhs_, "rhs"); CHKERRV(ierr);
  
  // initialize size of vector
  ierr = VecSetSizes(rhs_, PETSC_DECIDE, n); CHKERRV(ierr);
  
  // set sparsity type and other options
  ierr = VecSetFromOptions(rhs_);  CHKERRV(ierr);
  
  // initialize solution vector as copy of rhs
  ierr = VecDuplicate(rhs_, &solution_); CHKERRV(ierr);
  
  // create PETSc matrix object
  
  // PETSc MatCreateAIJ parameters
  int d_nz = 3;   // number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
  int o_nz = 0;   //  number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows)
  
  switch (mesh_->dimension())
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
  
  LOG(DEBUG) << "d="<<mesh_->dimension()
    <<", number of diagonal non-zeros: "<<d_nz<<", number of off-diagonal non-zeros: "<<o_nz;
  
  ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, n, n, 
                      d_nz, NULL, o_nz, NULL, &stiffnessMatrix_); CHKERRV(ierr);
  ierr = MatMPIAIJSetPreallocation(stiffnessMatrix_, d_nz, NULL, o_nz, NULL); CHKERRV(ierr);
  
  // dense matrix:
  if (false)
  {
    ierr = MatCreate(PETSC_COMM_WORLD, &stiffnessMatrix_);  CHKERRV(ierr);
    ierr = MatSetSizes(stiffnessMatrix_, PETSC_DECIDE,PETSC_DECIDE, n, n);  CHKERRV(ierr);
    ierr = MatSetFromOptions(stiffnessMatrix_); CHKERRV(ierr);
    ierr = MatSetUp(stiffnessMatrix_); CHKERRV(ierr);
  }
}

void FiniteElements::finalAssembly()
{
  PetscErrorCode ierr;
  // communicate portions to the right processors before using the matrix and vector in computations
  ierr = VecAssemblyBegin(rhs_); CHKERRV(ierr);
  ierr = VecAssemblyEnd(rhs_); CHKERRV(ierr);
  
  ierr = MatAssemblyBegin(stiffnessMatrix_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(stiffnessMatrix_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  
  LOG(DEBUG) << "finalAssembly";
}

Mat &FiniteElements::stiffnessMatrix()
{
  return stiffnessMatrix_;
}

Vec &FiniteElements::rightHandSide()
{
  return rhs_;
}

Vec &FiniteElements::solution()
{
  return solution_;
}

Mat &FiniteElements::discretizationMatrix()
{
  return discretizationMatrix_;
}

void FiniteElements::print()
{
  if (disablePrinting_)
    return;
  
  VLOG(1)<<"======================";
  int nRows, nColumns;
  MatGetSize(stiffnessMatrix_, &nRows, &nColumns);
  VLOG(1)<<"stiffnessMatrix ("<<nRows<<" x "<<nColumns<<") and rhs:";
  
  if (!disableMatrixPrinting_)
  {
    VLOG(1) << std::endl<<PetscUtility::getStringMatrixVector(stiffnessMatrix_, rhs_);
    VLOG(1) << "sparsity pattern: " << std::endl << PetscUtility::getStringSparsityPattern(stiffnessMatrix_);
  }
  
  MatInfo info;
  MatGetInfo(stiffnessMatrix_, MAT_LOCAL, &info);
  
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
  VecGetSize(rhs_, &nEntries);
  VLOG(1)<<"rhs ("<<nEntries<<" entries):";
  VLOG(1)<<PetscUtility::getStringVector(rhs_);
  VLOG(1)<<"======================";
  
  VecGetSize(solution_, &nEntries);
  VLOG(1)<<"solution ("<<nEntries<<" entries):";
  VLOG(1)<<PetscUtility::getStringVector(solution_);
  VLOG(1)<<"======================";
}

bool FiniteElements::discretizationMatrixInitialized()
{
  return discretizationMatrixInitialized_;
}

void FiniteElements::initializeDiscretizationMatrix()
{
  // determine problem size
  int nEntries;
  VecGetSize(rhs_, &nEntries);
  
  // create PETSc matrix object
  
  // PETSc MatCreateAIJ parameters
  int d_nz = 3;   // number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
  int o_nz = 0;   //  number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows)
  
  switch (mesh_->dimension())
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
                      d_nz, NULL, o_nz, NULL, &discretizationMatrix_); CHKERRV(ierr);
  ierr = MatMPIAIJSetPreallocation(discretizationMatrix_, d_nz, NULL, o_nz, NULL); CHKERRV(ierr);
  
  discretizationMatrixInitialized_ = true;
}

} // namespace Data
