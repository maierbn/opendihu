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
using FiniteElementsSolidMechanics = FiniteElements<
  BasisOnMeshType,
  Term,
  Equation::isSolidMechanics<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>;
 
template<typename BasisOnMeshType,typename Term>
FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
FiniteElements(DihuContext context) : Data<BasisOnMeshType>(context)
{
  LOG(TRACE) << "Data::FiniteElements constructor";
  PythonUtility::printDict(this->context_.getPythonConfig());
}

template<typename BasisOnMeshType,typename Term>
FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
~FiniteElements()
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
  return this->mesh_->nNodes() * 3;  // 3 components for displacements
}

template<typename BasisOnMeshType,typename Term>
void FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
createPetscObjects()
{
  dof_no_t n = this->mesh_->nDofs();
  
  LOG(DEBUG)<<"FiniteElementsSolidMechanics<BasisOnMeshType,Term>::createPetscObjects("<<n<<")";
  
  this->residual_ = this->mesh_->template createFieldVariable<3>("residual");
  this->externalVirtualEnergy_ = this->mesh_->template createFieldVariable<3>("externalVirtualEnergy");
  this->increment_ = this->mesh_->template createFieldVariable<3>("increment");
  this->displacements_ = this->mesh_->template createFieldVariable<3>("displacements");
  this->geometryReference_ = this->mesh_->template createFieldVariable<3>("geometryReference");
  
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
                      diagonalNonZeros, NULL, offdiagonalNonZeros, NULL, &this->tangentStiffnessMatrix_); CHKERRV(ierr);
    ierr = MatMPIAIJSetPreallocation(this->tangentStiffnessMatrix_, diagonalNonZeros, NULL, offdiagonalNonZeros, NULL); CHKERRV(ierr);
  }
  // allow additional non-zero entries in the stiffness matrix for UnstructuredDeformable mesh
  //MatSetOption(this->tangentStiffnessMatrix_, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
  
  // dense matrix
  if (false)
  {
    ierr = MatCreate(PETSC_COMM_WORLD, &this->tangentStiffnessMatrix_);  CHKERRV(ierr);
    ierr = MatSetSizes(this->tangentStiffnessMatrix_, PETSC_DECIDE,PETSC_DECIDE, n, n);  CHKERRV(ierr);
    ierr = MatSetFromOptions(this->tangentStiffnessMatrix_); CHKERRV(ierr);
    ierr = MatSetUp(this->tangentStiffnessMatrix_); CHKERRV(ierr);
  }
}

template<typename BasisOnMeshType,typename Term>
void FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
finalAssembly()
{
  PetscErrorCode ierr;
  // communicate portions to the right processors before using the matrix and vector in computations
  ierr = MatAssemblyBegin(this->tangentStiffnessMatrix_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  
  ierr = MatAssemblyEnd(this->tangentStiffnessMatrix_, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  
  LOG(DEBUG) << "finalAssembly";
}

template<typename BasisOnMeshType,typename Term>
Mat &FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
tangentStiffnessMatrix()
{
  return this->tangentStiffnessMatrix_;
}

template<typename BasisOnMeshType,typename Term>
FieldVariable::FieldVariable<BasisOnMeshType,3> &FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
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
FieldVariable::FieldVariable<BasisOnMeshType,3> &FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
displacements()
{
  return *this->displacements_;
}


template<typename BasisOnMeshType,typename Term>
FieldVariable::FieldVariable<BasisOnMeshType,3> &FiniteElementsSolidMechanics<BasisOnMeshType,Term>::
externalVirtualEnergy()
{
  return *this->externalVirtualEnergy_;
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
    externalVirtualEnergy_
  );
}
  

} // namespace Data
