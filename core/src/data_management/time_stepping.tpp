#include "data_management/time_stepping.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <numeric>
#include <memory>

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "control/dihu_context.h"
#include "utility/petsc_utility.h"

namespace Data
{

template<typename BasisOnMeshType,int nComponents>
TimeStepping<BasisOnMeshType,nComponents>::
TimeStepping(DihuContext context) : Data<BasisOnMeshType>(context)
{
}

template<typename BasisOnMeshType,int nComponents>
TimeStepping<BasisOnMeshType,nComponents>::
~TimeStepping()
{
  // free PETSc objects
  if (this->initialized_)
  {
    //PetscErrorCode ierr;
    //ierr = VecDestroy(&solution_); CHKERRV(ierr);
  }
}

template<typename BasisOnMeshType,int nComponents>
void TimeStepping<BasisOnMeshType,nComponents>::
createPetscObjects()
{
  LOG(DEBUG)<<"TimeStepping<BasisOnMeshType,nComponents>::createPetscObjects("<<nComponents<<")"<<std::endl;
  assert(this->mesh_);
  this->solution_ = this->mesh_->template createFieldVariable<nComponents>("solution");
  this->increment_ = this->mesh_->template createFieldVariable<nComponents>("increment");
}

template<typename BasisOnMeshType,int nComponents>
void TimeStepping<BasisOnMeshType,nComponents>::
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

template<typename BasisOnMeshType,int nComponents>
void TimeStepping<BasisOnMeshType,nComponents>::
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

template<typename BasisOnMeshType,int nComponents>
bool TimeStepping<BasisOnMeshType,nComponents>::
invLumMassMatrixInitialized()
{
  return this.invLumMassMatrixInitialized_;
}

template<typename BasisOnMeshType,int nComponents>
void TimeStepping<BasisOnMeshType,nComponents>::
initializeInvLumMassMatrix()
{
  // determine problem size
  int nEntries;
  VecGetSize(this->solution_->values(), &nEntries);
  
  // create PETSc matrix object
  
  // PETSc MatCreateAIJ parameters
  int diagonalNonZeros = 1;   // number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
  int offdiagonalNonZeros = 0;   //  number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows)
  
  getPetscMemoryParameters(diagonalNonZeros, offdiagonalNonZeros);
  
  PetscErrorCode ierr;
  ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, nEntries, nEntries,
                      diagonalNonZeros, NULL, offdiagonalNonZeros, NULL, &this->invLumMassMatrix_); CHKERRV(ierr);
                      ierr = MatMPIAIJSetPreallocation(this->massMatrix_, diagonalNonZeros, NULL, offdiagonalNonZeros, NULL); CHKERRV(ierr);
                      
                      this->invLumMassMatrixInitialized_ = true;
}

template<typename BasisOnMeshType,int nComponents>
Mat &TimeStepping<BasisOnMeshType,nComponents>::
invLumMassMatrix()
{
  return this->invLumMassMatrix_;
}

template<typename BasisOnMeshType,int nComponents>
Mat &TimeStepping<BasisOnMeshType,nComponents>::
systemMatrix()
{
  return this->systemMatrix_;
}

template<typename BasisOnMeshType,int nComponents>
FieldVariable::FieldVariable<BasisOnMeshType,nComponents> &TimeStepping<BasisOnMeshType,nComponents>::
solution()
{
  return *this->solution_;
}

template<typename BasisOnMeshType,int nComponents>
FieldVariable::FieldVariable<BasisOnMeshType,nComponents> &TimeStepping<BasisOnMeshType,nComponents>::
increment()
{
  return *this->increment_;
}

template<typename BasisOnMeshType,int nComponents>
dof_no_t TimeStepping<BasisOnMeshType,nComponents>::
nUnknowns()
{
  return this->mesh_->nNodes() * nComponents;
}

template<typename BasisOnMeshType,int nComponents>
constexpr int TimeStepping<BasisOnMeshType,nComponents>::
getNDofsPerNode()
{
  return nComponents;
}

template<typename BasisOnMeshType,int nComponents>
void TimeStepping<BasisOnMeshType,nComponents>::
print()
{
  if (!VLOG_IS_ON(4))
    return;

  VLOG(4)<<"======================";

  int nEntries;
  VecGetSize(this->increment_->values(), &nEntries);
  VLOG(4)<<"increment ("<<nEntries<<" entries):";
  VLOG(4)<<PetscUtility::getStringVector(this->increment_->values());
  VLOG(4)<<"======================";

  VecGetSize(this->solution_->values(), &nEntries);
  VLOG(4)<<"solution ("<<nEntries<<" entries):";
  VLOG(4)<<PetscUtility::getStringVector(this->solution_->values());
  VLOG(4)<<"======================";
}

template<typename BasisOnMeshType,int nComponents>
typename TimeStepping<BasisOnMeshType,nComponents>::OutputFieldVariables TimeStepping<BasisOnMeshType,nComponents>::
getOutputFieldVariables()
{
  return OutputFieldVariables(
    std::make_shared<FieldVariable::FieldVariable<BasisOnMeshType,3>>(this->mesh_->geometryField()),
    solution_
  );
}


} // namespace