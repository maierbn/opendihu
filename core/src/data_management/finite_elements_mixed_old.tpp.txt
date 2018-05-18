#include "data_management/finite_elements_mixed.h"

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

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
FiniteElements(DihuContext context) : Data<HighOrderBasisOnMeshType>(context)
{
  LOG(TRACE) << "Data::FiniteElements<Mixed> constructor";
  //PythonUtility::printDict(this->context_.getPythonConfig());
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
~FiniteElements()
{
  PetscErrorCode ierr;
  // free PETSc objects
  if (this->initialized_)
  {
    ierr = MatDestroy(&this->stiffnessMatrix_); CHKERRV(ierr);
  }
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
void FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
getPetscMemoryParameters(int &diagonalNonZeros, int &offdiagonalNonZeros)
{
  const int D = this->mesh_->dimension();
  const int nOverlaps = BasisOnMesh::BasisOnMeshBaseDim<1,typename HighOrderBasisOnMeshType::BasisFunction>::nDofsPerNode()*3;
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

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
void FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
createPetscObjects()
{
  dof_no_t n = this->mesh_->nDofs();
  
  LOG(DEBUG)<<"FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::createPetscObjects("<<n<<")";
  
  this->rhs_ = this->mesh_->template createFieldVariable<3>("rhs");
  this->solution_ = this->mesh_->template createFieldVariable<3>("solution");
  
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

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
void FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
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

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
Mat &FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
stiffnessMatrix()
{
  return this->stiffnessMatrix_;
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
FieldVariable::FieldVariable<HighOrderBasisOnMeshType,HighOrderBasisOnMeshType::dim()> &FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
rightHandSide()
{
  return *this->rhs_;
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
FieldVariable::FieldVariable<HighOrderBasisOnMeshType,HighOrderBasisOnMeshType::dim()> &FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
solution()
{
  return *this->solution_;
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
FieldVariable::FieldVariable<HighOrderBasisOnMeshType,HighOrderBasisOnMeshType::dim()> &FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
geometryReference()
{
  return *this->geometryReference_;
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
FieldVariable::FieldVariable<HighOrderBasisOnMeshType,HighOrderBasisOnMeshType::dim()> &FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
displacement()
{
  return *this->displacement_;
}


template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
FieldVariable::FieldVariable<LowOrderBasisOnMeshType,1> &FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
pressure()
{
  return *this->pressure_;
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
FieldVariable::FieldVariable<HighOrderBasisOnMeshType,1> &FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
f()
{
  return *this->f_;
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
Mat &FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
massMatrix()
{
  return this->massMatrix_;
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
Mat &FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
kuu()
{
  return kuu_;
} 

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
Mat &FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
kup()
{
  return kup_;
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
Mat &FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
kpp()
{
  return kpp_;
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
Mat &FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
kppInverse()
{
  return kppInverse_;
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
Mat &FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
kupTranspose()
{
  return kupTranspose_;
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
Mat &FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
tempKupMatrix()
{
  return tempKupMatrix_;
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
Mat &FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
schurComplement()
{
  return schurComplement_;
}
    
template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
void FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
getLocalVectors(Vec &fu, Vec &fp, Vec &tempKppFp, Vec &tempKupKppFp)
{
  fu = fu_;
  fp = fp_;
  tempKppFp = tempKppFp_;
  tempKupKppFp = tempKupKppFp_;
}
    
template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
void FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
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

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
bool FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
massMatrixInitialized()
{
  return this->massMatrixInitialized_;
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
void FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
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

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
typename FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::OutputFieldVariables FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
getOutputFieldVariables()
{
  return OutputFieldVariables(
    std::make_shared<FieldVariable::FieldVariable<HighOrderBasisOnMeshType,3>>(this->mesh_->geometryField()),
    solution_,
    rhs_
  );
}
  
template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
void FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
setMesh(std::shared_ptr<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>> mixedMesh)
{
  mixedMesh_ = mixedMesh;
  
  // store high order mesh as mesh_
  this->mesh_ = mixedMesh_->highOrderBasisOnMesh();
} 

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
void FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
initialize()
{
  Data<HighOrderBasisOnMeshType>::initialize();
 
  LOG(DEBUG) << "mesh has geometry field: " << this->mesh_->hasGeometryField();
  initializeFieldVariables();
  initializeMatrices();
  
  // set up diffusion tensor if there is any
  DiffusionTensor<HighOrderBasisOnMeshType::dim()>::initialize(this->context_.getPythonConfig());
}
  
template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
void FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
initializeFieldVariables()
{
  // generate geometryReference variable as copy of geometry field
  assert(this->mesh_->hasGeometryField());
  geometryReference_ = std::make_shared<FieldVariable::FieldVariable<HighOrderBasisOnMeshType,3>>(this->mesh_->geometryField(), "geometryReference");
  //geometryReference_->initializeFromFieldVariable(this->mesh_->geometryField(), "geometryReference", {"x","y","z"});
  //geometryReference_->setValues(this->mesh_->geometryField());
  
  std::vector<std::string> components({"x","y","z"});
  displacement_ = std::make_shared<FieldVariable::FieldVariable<HighOrderBasisOnMeshType,3>>(this->mesh_, "displacement", components);
  //displacement_->initializeFromFieldVariable(this->mesh_->geometryField(), "displacement", {"x","y","z"});
  
  std::vector<std::string> unnamedSingleComponent({"0"});
  pressure_ = std::make_shared<FieldVariable::FieldVariable<LowOrderBasisOnMeshType,1>>(this->mixedMesh_->lowOrderBasisOnMesh(), "pressure", unnamedSingleComponent);
  
  f_ = std::make_shared<FieldVariable::FieldVariable<HighOrderBasisOnMeshType,1>>(this->mixedMesh_->highOrderBasisOnMesh(), "f", unnamedSingleComponent);
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
void FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
initializeMatrices()
{
  PetscErrorCode ierr;
  
  // create PETSc matrix objects
  ierr = MatCreate(PETSC_COMM_WORLD, &this->kuu_);  CHKERRV(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD, &this->kup_);  CHKERRV(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD, &this->kpp_);  CHKERRV(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD, &this->kppInverse_);  CHKERRV(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD, &this->kupTranspose_);  CHKERRV(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD, &this->schurComplement_);  CHKERRV(ierr);
  
  const int nu = HighOrderBasisOnMeshType::nDofsPerElement();
  const int np = LowOrderBasisOnMeshType::nDofsPerElement();
  ierr = MatSetSizes(this->kuu_, nu, nu, nu, nu);  CHKERRV(ierr);
  ierr = MatSetUp(this->kuu_); CHKERRV(ierr);
  
  ierr = MatSetSizes(this->kup_, nu, np, nu, np);  CHKERRV(ierr);
  ierr = MatSetUp(this->kup_); CHKERRV(ierr);
  
  ierr = MatSetSizes(this->kpp_, np, np, np, np);  CHKERRV(ierr);
  ierr = MatSetUp(this->kpp_); CHKERRV(ierr);
  
  ierr = MatSetSizes(this->kpp_, np, np, np, np);  CHKERRV(ierr);
  ierr = MatSetUp(this->kppInverse_); CHKERRV(ierr);
  
  ierr = MatSetSizes(this->kupTranspose_, np, nu, np, nu);  CHKERRV(ierr);
  ierr = MatSetUp(this->kupTranspose_); CHKERRV(ierr);
  
  ierr = MatSetSizes(this->tempKupMatrix_, nu, np, nu, np);  CHKERRV(ierr);
  ierr = MatSetUp(this->tempKupMatrix_); CHKERRV(ierr);
  
  ierr = MatSetSizes(this->schurComplement_, nu, nu, nu, nu);  CHKERRV(ierr);
  ierr = MatSetUp(this->schurComplement_); CHKERRV(ierr);
  
  // create PETSc vector objects
  ierr = VecCreate(PETSC_COMM_WORLD, &this->fu_);  CHKERRV(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD, &this->fp_);  CHKERRV(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD, &this->tempKppFp_);  CHKERRV(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD, &this->tempKupKppFp_);  CHKERRV(ierr);
  
  ierr = VecSetSizes(this->fu_, nu, nu); CHKERRV(ierr);
  ierr = VecSetSizes(this->fp_, np, np); CHKERRV(ierr);
  ierr = VecSetSizes(this->tempKppFp_, np, np); CHKERRV(ierr);
  ierr = VecSetSizes(this->tempKupKppFp_, nu, nu); CHKERRV(ierr);
  VecSetType(this->fu_, VECSEQ);
  VecSetType(this->fp_, VECSEQ);
  VecSetType(this->tempKppFp_, VECSEQ);
  VecSetType(this->tempKupKppFp_, VECSEQ);
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType,typename Term>
std::shared_ptr<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>> FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>,Term>::
mixedMesh()
{
  return mixedMesh_;
}



} // namespace Data
