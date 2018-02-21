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
#include "basis_on_mesh/05_basis_on_mesh.h"
#include "mesh/unstructured_deformable.h"
#include "basis_function/hermite.h"

namespace Data
{

template<typename BasisOnMeshType>
FiniteElements<BasisOnMeshType>::
FiniteElements(const DihuContext &context) : Data<BasisOnMeshType>(context)
{
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

/*
template<typename BasisOnMeshType>
void FiniteElements<BasisOnMeshType>::
getPetscMemoryParameters(int &diagonalNonZeros, int &offdiagonalNonZeros)
{
  const int nOverlaps = 3;
  switch (this->mesh_->dimension())
  {
  case 1:
    diagonalNonZeros = nOverlaps;
    offdiagonalNonZeros = nOverlaps;
    break;
  case 2:
    diagonalNonZeros = pow(nOverlaps, 2);
    offdiagonalNonZeros = diagonalNonZeros;
    break;
  case 3:
    diagonalNonZeros = pow(nOverlaps, 3);
    offdiagonalNonZeros = diagonalNonZeros;
    break;
  };
}*/

// for UnstructuredDeformable and Hermite
//template<int D>
//void FiniteElements<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>, BasisFunction::Hermite>>::


//template<int D, typename BasisFunctionType>
//void FiniteElements<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>::

template<typename BasisOnMeshType>
void FiniteElements<BasisOnMeshType>::
getPetscMemoryParameters(int &diagonalNonZeros, int &offdiagonalNonZeros)
{
  const int D = this->mesh_->dimension();
  const int nOverlaps = BasisOnMesh::BasisOnMeshBaseDim<1,typename BasisOnMeshType::BasisFunction>::nDofsPerNode()*3;
  switch (D)
  {
  case 1:
    diagonalNonZeros = nOverlaps;
    offdiagonalNonZeros = nOverlaps;
    break;
  case 2:
    diagonalNonZeros = pow(nOverlaps, 2);
    offdiagonalNonZeros = diagonalNonZeros;
    break;
  case 3:
    diagonalNonZeros = pow(nOverlaps, 3);
    offdiagonalNonZeros = diagonalNonZeros;
    break;
  };
}

template<typename BasisOnMeshType>
void FiniteElements<BasisOnMeshType>::
createPetscObjects()
{
  dof_no_t n = this->mesh_->nDofs();
  
  LOG(DEBUG)<<"FiniteElements<BasisOnMeshType>::createPetscObjects("<<n<<")";
  
  this->rhs_ = this->mesh_->createFieldVariable("rhs");
  this->solution_ = this->mesh_->createFieldVariable("solution");
  
  PetscErrorCode ierr;
  // create PETSc matrix object
  
  // PETSc MatCreateAIJ parameters
  int diagonalNonZeros = 3;   // number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
  int offdiagonalNonZeros = 0;   //  number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows)
  
  getPetscMemoryParameters(diagonalNonZeros, offdiagonalNonZeros);
  
  LOG(DEBUG) << "d="<<this->mesh_->dimension()
    <<", number of diagonal non-zeros: "<<diagonalNonZeros<<", number of off-diagonal non-zeros: "<<offdiagonalNonZeros;
  
  // sparse matrix
  if (true)
  {
    ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, n, n, 
                      diagonalNonZeros, NULL, offdiagonalNonZeros, NULL, &this->stiffnessMatrix_); CHKERRV(ierr);
    ierr = MatMPIAIJSetPreallocation(this->stiffnessMatrix_, diagonalNonZeros, NULL, offdiagonalNonZeros, NULL); CHKERRV(ierr);
  }
  // allow additional non-zero entries in the stiffness matrix for UnstructuredDeformable mesh
  //MatSetOption(this->stiffnessMatrix_, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
  
  // dense matrix
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
massMatrix()
{
  return this->massMatrix_;
}

template<typename BasisOnMeshType>
void FiniteElements<BasisOnMeshType>::
print()
{
  if (!VLOG_IS_ON(4))
    return;
  
  VLOG(4)<<"======================";
  int nRows, nColumns;
  MatGetSize(this->stiffnessMatrix_, &nRows, &nColumns);
  VLOG(4)<<"stiffnessMatrix ("<<nRows<<" x "<<nColumns<<") and rhs:";
  
  if (!this->disableMatrixPrinting_)
  {
    VLOG(4) << std::endl<<PetscUtility::getStringMatrixVector(this->stiffnessMatrix_, this->rhs_->values());
    VLOG(4) << "sparsity pattern: " << std::endl << PetscUtility::getStringSparsityPattern(this->stiffnessMatrix_);
  }
  
  MatInfo info;
  MatGetInfo(this->stiffnessMatrix_, MAT_LOCAL, &info);
  
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
  VecGetSize(this->rhs_->values(), &nEntries);
  VLOG(4)<<"rhs ("<<nEntries<<" entries):";
  VLOG(4)<<PetscUtility::getStringVector(this->rhs_->values());
  VLOG(4)<<"======================";
  
  VecGetSize(this->solution_->values(), &nEntries);
  VLOG(4)<<"solution ("<<nEntries<<" entries):";
  VLOG(4)<<PetscUtility::getStringVector(this->solution_->values());
  VLOG(4)<<"======================";
}

template<typename BasisOnMeshType>
bool FiniteElements<BasisOnMeshType>::
massMatrixInitialized()
{
  return this->massMatrixInitialized_;
}

template<typename BasisOnMeshType>
void FiniteElements<BasisOnMeshType>::
initializeMassMatrix()
{
  // determine problem size
  int nEntries;
  VecGetSize(this->rhs_->values(), &nEntries);
  
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

template<typename BasisOnMeshType>
std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> FiniteElements<BasisOnMeshType>::
fieldVariables()
{
  std::vector<std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType>>> result;
  result.push_back(std::make_shared<FieldVariable::FieldVariable<BasisOnMeshType>>(this->mesh_->geometryField()));
  result.push_back(solution_);
  result.push_back(rhs_);
  this->mesh_->addNonGeometryFieldVariables(result);   // add all further field variables that were e.g. present in an input file
  
  return result;
}
  

} // namespace Data
