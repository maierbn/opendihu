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

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
FiniteElements(DihuContext context) : Data<BasisOnMeshType>(context)
{
  LOG(TRACE) << "Data::FiniteElements constructor";
  //PythonUtility::printDict(this->context_.getPythonConfig());
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
~FiniteElements()
{
  PetscErrorCode ierr;
  // free PETSc objects
  if (this->initialized_)
  {
    ierr = MatDestroy(&this->stiffnessMatrix_); CHKERRV(ierr);
  }
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
void FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
initialize()
{
  Data<BasisOnMeshType>::initialize();

  // set up diffusion tensor if there is any
  DiffusionTensor<BasisOnMeshType::dim()>::initialize(this->context_.getPythonConfig());
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
void FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
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

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
void FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
createPetscObjects()
{
  FieldVariable::FieldVariable<BasisOnMeshType,3> geometry = this->mesh_->geometryField();
 
  if (VLOG_IS_ON(1))
  {
    //VLOG(1) << PetscUtility::getStringVector(geometry.values());
  }
  
  
  dof_no_t n = this->mesh_->nDofs();

  LOG(DEBUG)<<"FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::createPetscObjects("<<n<<")";

  this->rhs_ = this->mesh_->template createFieldVariable<1>("rhs");
  this->solution_ = this->mesh_->template createFieldVariable<1>("solution");

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
  
  createPetscObjects_systemMatrix();
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
void FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
createPetscObjects_systemMatrix()
{
  dof_no_t n = this->mesh_->nDofs();
  
  LOG(DEBUG)<<"TimeStepping<BasisOnMeshType,nComponents>::createPetscObjects_systemMatrix("<<n<<")";
  
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
                        diagonalNonZeros, NULL, offdiagonalNonZeros, NULL, &this->systemMatrix_); CHKERRV(ierr);
                        ierr = MatMPIAIJSetPreallocation(this->systemMatrix_, diagonalNonZeros, NULL, offdiagonalNonZeros, NULL); CHKERRV(ierr);
  }
  // allow additional non-zero entries in the stiffness matrix for UnstructuredDeformable mesh
  //MatSetOption(this->stiffnessMatrix_, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
  
  // dense matrix
  if (false)
  {
    ierr = MatCreate(PETSC_COMM_WORLD, &this->systemMatrix_);  CHKERRV(ierr);
    ierr = MatSetSizes(this->systemMatrix_, PETSC_DECIDE,PETSC_DECIDE, n, n);  CHKERRV(ierr);
    ierr = MatSetFromOptions(this->systemMatrix_); CHKERRV(ierr);
    ierr = MatSetUp(this->systemMatrix_); CHKERRV(ierr);
  }
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
void FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
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

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
Mat &FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
stiffnessMatrix()
{
  return this->stiffnessMatrix_;
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
Mat &FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
systemMatrix()
{
  return this->systemMatrix_;
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
FieldVariable::FieldVariable<BasisOnMeshType,1> &FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
rightHandSide()
{
  return *this->rhs_;
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
FieldVariable::FieldVariable<BasisOnMeshType,1> &FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
solution()
{
  return *this->solution_;
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
Mat &FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
massMatrix()
{
  return this->massMatrix_;
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
Mat &FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
invLumMassMatrix()
{
  return this->invLumMassMatrix_;
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
void FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
print()
{
  if (!VLOG_IS_ON(4))
    return;

  VLOG(4)<<"======================";
  int nRows, nColumns;
  MatGetSize(this->stiffnessMatrix_, &nRows, &nColumns);
  VLOG(4)<<"stiffnessMatrix ("<<nRows<<" x "<<nColumns<<") and rhs:";

  VLOG(4) << std::endl<<PetscUtility::getStringMatrixVector(this->stiffnessMatrix_, this->rhs_->values());
  VLOG(4) << "sparsity pattern: " << std::endl << PetscUtility::getStringSparsityPattern(this->stiffnessMatrix_);

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

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
bool FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
massMatrixInitialized()
{
  return this->massMatrixInitialized_;
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
bool FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
invLumMassMatrixInitialized()
{
  return this->invLumMassMatrixInitialized_;
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
void FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
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

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
void FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
initializeInvLumMassMatrix()
{
  const int D = BasisOnMeshType::dim();
  LOG(TRACE)<<"initializeInvLumMassMatrix" << D << "D";
  
  // determine problem size
  int nEntries;
  VecGetSize(this->solution_->values(), &nEntries);
  LOG(INFO)<<"nEntries " << nEntries;
  LOG(TRACE)<<"checkpoint";
  
  // create PETSc matrix object
  
  // PETSc MatCreateAIJ parameters
  int diagonalNonZeros = 3;   // number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
  int offdiagonalNonZeros = 0;   //  number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows)
  
  getPetscMemoryParameters(diagonalNonZeros, offdiagonalNonZeros);
  
  PetscErrorCode ierr;
  ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, nEntries, nEntries,
                      diagonalNonZeros, NULL, offdiagonalNonZeros, NULL, &this->invLumMassMatrix_); CHKERRV(ierr);
  ierr = MatMPIAIJSetPreallocation(this->invLumMassMatrix_, diagonalNonZeros, NULL, offdiagonalNonZeros, NULL); CHKERRV(ierr);
                      
  this->invLumMassMatrixInitialized_ = true;
  
  LOG(TRACE)<<"finished";
}

template<typename BasisOnMeshType,typename Term,typename DummyForTraits,typename DummyForTraits2>
typename FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::OutputFieldVariables FiniteElements<BasisOnMeshType,Term,DummyForTraits,DummyForTraits2>::
getOutputFieldVariables()
{
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>> geometryField
    = std::make_shared<FieldVariable::FieldVariable<BasisOnMeshType,3>>(this->mesh_->geometryField());
  /*
  std::shared_ptr<FieldVariable::FieldVariable<BasisOnMeshType,3>> generalField;

  generalField = std::static_pointer_cast<FieldVariable::FieldVariable<BasisOnMeshType,3>>(this->mesh_->fieldVariable("general"));
  if (!generalField)
   generalField = geometryField;
  */
  return OutputFieldVariables(
    geometryField,
    solution_,
    rhs_
  );
}


} // namespace Data
