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

FiniteElements::FiniteElements(DihuContext &context) : Data(context)
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
  element_idx_t n = mesh_->nDegreesOfFreedom();
  
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
    d_nz = 24;
    o_nz = 24;
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
  
  LOG(INFO)<<"======================";
  int nRows, nColumns;
  MatGetSize(stiffnessMatrix_, &nRows, &nColumns);
  
  LOG(INFO)<<"stiffnessMatrix ("<<nRows<<" x "<<nColumns<<") and rhs:";
  
  if (!disableMatrixPrinting_)
  {
    std::vector<double> matrixValues; 
    std::vector<double> vectorValues;
    
    PetscUtility::getMatrixEntries(stiffnessMatrix_, matrixValues);
    PetscUtility::getVectorEntries(rhs_, vectorValues);
    
    std::stringstream s;
    s<<"    ";
    for (int j=0; j<nColumns; j++)
    {
      s<<std::setw(5)<<std::setfill('_')<<j;
    }
    s<<std::string(5,'_')<<" | rhs";
    s<<std::endl;
    for (int i=0; i<nRows; i++)
    {
      s<<std::setw(3)<<std::setfill(' ')<<i<<"| ";
      for (int j=0; j<nColumns; j++)
      {
        if(matrixValues[i*nRows + j] == 0.0)
          s<<std::string(5, ' ');
        else
          s<<std::setw(4)<<std::setfill(' ')<<matrixValues[i*nRows + j]<<" ";
      }
      s<<std::string(5, ' ')<<"| "<<vectorValues[i];
      s<<std::endl;
    }
    s<<std::endl;
    LOG(INFO) << std::endl<<s.str();
    
    
    s.str("");
    s<<" ";
    for (int j=0; j<nColumns; j++)
    {
      if (j%10 == 0)
        s<<"|";
      else if (j%2 == 0)
        s<<".";
      else 
        s<<" ";
    }
    s<<std::endl;
    for (int i=0; i<nRows; i++)
    {
      s<<" ";
      for (int j=0; j<nColumns; j++)
      {
        if(matrixValues[i*nRows + j] == 0.0)
          s<<" ";
        else
          s<<"*";
      }
      s<<std::endl;
    }
    s<<std::endl;
    LOG(INFO) << "sparsity pattern: " << std::endl<<s.str();
  }
  
  MatInfo info;
  MatGetInfo(stiffnessMatrix_, MAT_LOCAL, &info);
  
  LOG(INFO)<<"Matrix info: "<<std::endl
    <<"block_size: "<<info.block_size<<std::endl
    <<"number of nonzeros: allocated: "<<info.nz_allocated<<", used: "<<info.nz_used<<", unneeded: "<<info.nz_unneeded<<std::endl
    <<"memory allocated: "<<info.memory<<std::endl
    <<"number of matrix assemblies called: "<<info.assemblies<<std::endl
    <<"number of mallocs during MatSetValues(): "<<info.mallocs<<std::endl
    <<"fill ratio for LU/ILU: given: "<<info.fill_ratio_given<<", needed: "<<info.fill_ratio_needed<<std::endl 
    <<"number of mallocs during factorization: "<<info.factor_mallocs<<std::endl;
    
    
  LOG(INFO)<<"======================";
  
  int nEntries;
  VecGetSize(rhs_, &nEntries);
  
  LOG(INFO)<<"rhs ("<<nEntries<<" entries):";
  
  std::vector<double> vectorValues;
  PetscUtility::getVectorEntries(rhs_, vectorValues);
  
  std::stringstream s;
  for (int i=0; i<nEntries; i++)
  {
    s<<vectorValues[i]<<" ";
  }
  LOG(INFO)<<s.str();
  LOG(INFO)<<"======================";
  
  VecGetSize(solution_, &nEntries);
  
  LOG(INFO)<<"solution ("<<nEntries<<" entries):";
    
  PetscUtility::getVectorEntries(solution_, vectorValues);
  s.str("");
  for (int i=0; i<nEntries; i++)
  {
    s<<vectorValues[i]<<" ";
  }
  LOG(INFO)<<s.str();
  LOG(INFO)<<"======================";
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
    d_nz = 24;
    o_nz = 24;
    break;
  };
  
  PetscErrorCode ierr;
  ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, nEntries, nEntries, 
                      d_nz, NULL, o_nz, NULL, &discretizationMatrix_); CHKERRV(ierr);
  ierr = MatMPIAIJSetPreallocation(discretizationMatrix_, d_nz, NULL, o_nz, NULL); CHKERRV(ierr);
  
  discretizationMatrixInitialized_ = true;
}

} // namespace Data
