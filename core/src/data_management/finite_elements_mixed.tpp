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

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>>::
FiniteElements(const DihuContext &context) : Data<HighOrderBasisOnMeshType>(context)
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

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>>::
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
template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
void FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>>::
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
//void FiniteElements<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>, BasisFunction::Hermite>>::


//template<int D, typename BasisFunctionType>
//void FiniteElements<BasisOnMesh::BasisOnMesh<Mesh::UnstructuredDeformable<D>,BasisFunctionType>>::

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
void FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>>::
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

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
void FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>>::
createPetscObjects()
{
  dof_no_t n = this->mesh_->nDofs();
  
  LOG(DEBUG)<<"FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>>::createPetscObjects("<<n<<")";
  
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

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
void FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>>::
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

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
Mat &FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>>::
stiffnessMatrix()
{
  return this->stiffnessMatrix_;
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
FieldVariable::FieldVariable<HighOrderBasisOnMeshType> &FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>>::
rightHandSide()
{
  return *this->rhs_;
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
FieldVariable::FieldVariable<HighOrderBasisOnMeshType> &FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>>::
solution()
{
  return *this->solution_;
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
FieldVariable::FieldVariable<HighOrderBasisOnMeshType> &FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>>::
geometryReference()
{
  return *this->geometryReference_;
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
FieldVariable::FieldVariable<HighOrderBasisOnMeshType> &FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>>::
displacement()
{
  return *this->displacement_;
}


template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
FieldVariable::FieldVariable<LowOrderBasisOnMeshType> &FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>>::
pressure()
{
  return *this->pressure_;
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
Mat &FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>>::
discretizationMatrix()
{
  return this->discretizationMatrix_;
}

//! return the element stiffness matrix kuu
template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
Mat &FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>>::
kuu()
{
  return kuu_;
} 

//! return the element stiffness matrix kup
template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
Mat &FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>>::
kup()
{
  return kup_;
}

//! return the element stiffness matrix kpp
template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
Mat &FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>>::
kpp()
{
  return kpp_;
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
void FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>>::
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

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
bool FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>>::
discretizationMatrixInitialized()
{
  return this->discretizationMatrixInitialized_;
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
void FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>>::
initializeDiscretizationMatrix()
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
                      diagonalNonZeros, NULL, offdiagonalNonZeros, NULL, &this->discretizationMatrix_); CHKERRV(ierr);
  ierr = MatMPIAIJSetPreallocation(this->discretizationMatrix_, diagonalNonZeros, NULL, offdiagonalNonZeros, NULL); CHKERRV(ierr);
  
  this->discretizationMatrixInitialized_ = true;
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
std::vector<std::shared_ptr<FieldVariable::FieldVariable<HighOrderBasisOnMeshType>>> FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>>::
fieldVariables()
{
  std::vector<std::shared_ptr<FieldVariable::FieldVariable<HighOrderBasisOnMeshType>>> result;
  result.push_back(std::make_shared<FieldVariable::FieldVariable<HighOrderBasisOnMeshType>>(this->mesh_->geometryField()));
  result.push_back(solution_);
  result.push_back(rhs_);
  this->mesh_->addNonGeometryFieldVariables(result);   // add all further field variables that were e.g. present in an input file
  
  return result;
}
  
template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
void FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>>::
setMesh(std::shared_ptr<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>> mixedMesh)
{
  mixedMesh_ = mixedMesh;
  
  // store high order mesh as mesh_
  this->mesh_ = mixedMesh_->highOrderBasisOnMesh();
} 

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
void FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>>::
initialize()
{
  Data<HighOrderBasisOnMeshType>::initialize();
 
  LOG(DEBUG) << "mesh has geometry field: " << this->mesh_->hasGeometryField();
  initializeFieldVariables();
  initializeMatrices();
}
  
template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
void FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>>::
initializeFieldVariables()
{
  // generate geometryReference variable as copy of geometry field
  assert(this->mesh_->hasGeometryField());
  geometryReference_ = std::make_shared<FieldVariable::FieldVariable<HighOrderBasisOnMeshType>>(this->mesh_->geometryField(), "geometryReference");
  //geometryReference_->initializeFromFieldVariable(this->mesh_->geometryField(), "geometryReference", {"x","y","z"});
  //geometryReference_->setValues(this->mesh_->geometryField());
  
  std::vector<std::string> components({"x","y","z"});
  displacement_ = std::make_shared<FieldVariable::FieldVariable<HighOrderBasisOnMeshType>>(this->mesh_, "displacement", components);
  //displacement_->initializeFromFieldVariable(this->mesh_->geometryField(), "displacement", {"x","y","z"});
  
  pressure_ = std::make_shared<FieldVariable::FieldVariable<LowOrderBasisOnMeshType>>();
  
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
void FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>>::
initializeMatrices()
{
  PetscErrorCode ierr;
  
  // create PETSc matrix object
  ierr = MatCreate(PETSC_COMM_WORLD, &this->kuu_);  CHKERRV(ierr);
  
  const int nu = HighOrderBasisOnMeshType::nDofsPerElement();
  const int np = LowOrderBasisOnMeshType::nDofsPerElement();
  ierr = MatSetSizes(this->kuu_, nu, nu, nu, nu);  CHKERRV(ierr);
  ierr = MatSetUp(this->kuu_); CHKERRV(ierr);
  
  ierr = MatSetSizes(this->kup_, nu, np, nu, np);  CHKERRV(ierr);
  ierr = MatSetUp(this->kup_); CHKERRV(ierr);
  
  ierr = MatSetSizes(this->kpp_, np, np, np, np);  CHKERRV(ierr);
  ierr = MatSetUp(this->kpp_); CHKERRV(ierr);
}

template<typename LowOrderBasisOnMeshType,typename HighOrderBasisOnMeshType>
std::shared_ptr<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>> FiniteElements<BasisOnMesh::Mixed<LowOrderBasisOnMeshType,HighOrderBasisOnMeshType>>::
mixedMesh()
{
  return mixedMesh_;
}



} // namespace Data
