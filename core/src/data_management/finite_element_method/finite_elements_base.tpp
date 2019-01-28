#include "data_management/finite_element_method/finite_elements_base.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <memory>

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "control/dihu_context.h"
#include "utility/petsc_utility.h"
#include "function_space/function_space.h"
#include "mesh/unstructured_deformable.h"
#include "basis_function/hermite.h"
#include "partition/partitioned_petsc_mat/partitioned_petsc_mat.h"

namespace Data
{

template<typename FunctionSpaceType>
FiniteElementsBase<FunctionSpaceType>::
FiniteElementsBase(DihuContext context) : Data<FunctionSpaceType>(context)
{
}

template<typename FunctionSpaceType>
FiniteElementsBase<FunctionSpaceType>::
~FiniteElementsBase()
{
}

template<typename FunctionSpaceType>
void FiniteElementsBase<FunctionSpaceType>::
initialize()
{
  LOG(DEBUG) << "Data::FiniteElementsBase::initialize";
  Data<FunctionSpaceType>::initialize();
}
template<typename FunctionSpaceType>
void FiniteElementsBase<FunctionSpaceType>::
reset()
{
  LOG(DEBUG) << "Data::FiniteElementsBase::reset";
  // set initalize_ to false
  Data<FunctionSpaceType>::reset();

  // deallocate Petsc matrices
  this->stiffnessMatrix_ = nullptr;
}

template<typename FunctionSpaceType>
void FiniteElementsBase<FunctionSpaceType>::
getPetscMemoryParameters(int &diagonalNonZeros, int &offdiagonalNonZeros)
{
  const int D = this->functionSpace_->dimension();
  const int nDofsPerNode = FunctionSpace::FunctionSpaceBaseDim<1,typename FunctionSpaceType::BasisFunction>::nDofsPerNode();
  const int nDofsPerBasis = FunctionSpace::FunctionSpaceBaseDim<1,typename FunctionSpaceType::BasisFunction>::nDofsPerElement();
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

template<typename FunctionSpaceType>
void FiniteElementsBase<FunctionSpaceType>::
createPetscObjects()
{
  LOG(TRACE) << "FiniteElements::createPetscObjects";

  // get the partitioning from the function space
  std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> meshPartition = this->functionSpace_->meshPartition();
  
  // create field variables on local partition
  this->rhs_ = this->functionSpace_->template createFieldVariable<1>("rhs");
  this->solution_ = this->functionSpace_->template createFieldVariable<1>("solution");

  // create PETSc matrix object

  // PETSc MatCreateAIJ parameters
  int diagonalNonZeros = 3;   // number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
  int offdiagonalNonZeros = 0;   //  number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows)

  getPetscMemoryParameters(diagonalNonZeros, offdiagonalNonZeros);

  LOG(DEBUG) << "d=" << this->functionSpace_->dimension()
    << ", number of diagonal non-zeros: " << diagonalNonZeros << ", number of off-diagonal non-zeros: " <<offdiagonalNonZeros;

  int nComponents = 1;
  this->stiffnessMatrix_ = std::make_shared<PartitionedPetscMat<FunctionSpaceType>>(meshPartition, nComponents, diagonalNonZeros, offdiagonalNonZeros, "stiffnessMatrix");

}

template<typename FunctionSpaceType>
std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> FiniteElementsBase<FunctionSpaceType>::
stiffnessMatrix()
{
  return this->stiffnessMatrix_;
}

template<typename FunctionSpaceType>
std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> FiniteElementsBase<FunctionSpaceType>::
massMatrix()
{
  return this->massMatrix_;
}

template<typename FunctionSpaceType>
std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> FiniteElementsBase<FunctionSpaceType>::
inverseLumpedMassMatrix()
{
  return this->inverseLumpedMassMatrix_;
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> FiniteElementsBase<FunctionSpaceType>::
rightHandSide()
{
  return this->rhs_;
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> FiniteElementsBase<FunctionSpaceType>::
solution()
{
  return this->solution_;
}

template<typename FunctionSpaceType>
typename FiniteElementsBase<FunctionSpaceType>::TransferableSolutionDataType FiniteElementsBase<FunctionSpaceType>::
getSolutionForTransfer()
{
  LOG(INFO) << "retutrn FiniteElementsBase solution";
  return this->solution_;
}

template<typename FunctionSpaceType>
void FiniteElementsBase<FunctionSpaceType>::
print()
{
  if (!VLOG_IS_ON(4))
    return;

  VLOG(4) << "======================";
  VLOG(4) << *this->stiffnessMatrix_;
  VLOG(4) << *this->rhs_;
  VLOG(4) << *this->solution_;

  if (this->massMatrix_)
    VLOG(4) << *this->massMatrix_;

  if (this->inverseLumpedMassMatrix_)
    VLOG(4) << *this->inverseLumpedMassMatrix_;

  VLOG(4) << this->functionSpace_->geometryField();
  
  MatInfo info;
  MatGetInfo(this->stiffnessMatrix_->valuesGlobal(), MAT_LOCAL, &info);

  VLOG(4) << "stiffnessMatrix info: " << std::endl
    << "block_size: " <<info.block_size << std::endl
    << "number of nonzeros: allocated: " <<info.nz_allocated<< ", used: " <<info.nz_used<< ", unneeded: " <<info.nz_unneeded<< std::endl
    << "memory allocated: " <<info.memory<< std::endl
    << "number of matrix assemblies called: " <<info.assemblies << std::endl
    << "number of mallocs during MatSetValues(): " <<info.mallocs << std::endl
    << "fill ratio for LU/ILU: given: " <<info.fill_ratio_given<< ", needed: " <<info.fill_ratio_needed<< std::endl
    << "number of mallocs during factorization: " <<info.factor_mallocs << std::endl;

  VLOG(4) << "======================";
}

template<typename FunctionSpaceType>
void FiniteElementsBase<FunctionSpaceType>::
initializeMassMatrix()
{
  // if the massMatrix is already initialized do not initialize again
  if (this->massMatrix_)
    return;

  // create PETSc matrix object

  // PETSc MatCreateAIJ parameters
  int diagonalNonZeros = 3;   // number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
  int offdiagonalNonZeros = 0;   //  number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows)

  getPetscMemoryParameters(diagonalNonZeros, offdiagonalNonZeros);

  std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> partition = this->functionSpace_->meshPartition();
  const int nComponents = 1;
  this->massMatrix_ = std::make_shared<PartitionedPetscMat<FunctionSpaceType>>(partition, nComponents, diagonalNonZeros, offdiagonalNonZeros, "massMatrix");
}

template<typename FunctionSpaceType>
void FiniteElementsBase<FunctionSpaceType>::
initializeInverseLumpedMassMatrix()
{
  // if the inverseLumpedMassMatrix_ is already initialized do not initialize again
  if (this->inverseLumpedMassMatrix_)
    return;

  // create PETSc matrix object

  // PETSc MatCreateAIJ parameters
  int diagonalNonZeros = 3;   // number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
  int offdiagonalNonZeros = 0;   //  number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows)

  getPetscMemoryParameters(diagonalNonZeros, offdiagonalNonZeros);

  std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> partition = this->functionSpace_->meshPartition();
  const int nComponents = 1;
  this->inverseLumpedMassMatrix_ = std::make_shared<PartitionedPetscMat<FunctionSpaceType>>(partition, nComponents, diagonalNonZeros, offdiagonalNonZeros, "inverseLumpedMassMatrix");
}

template<typename FunctionSpaceType>
typename FiniteElementsBase<FunctionSpaceType>::OutputFieldVariables FiniteElementsBase<FunctionSpaceType>::
getOutputFieldVariables()
{
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> geometryField
    = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,3>>(this->functionSpace_->geometryField());
  /*
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> generalField;

  generalField = std::static_pointer_cast<FieldVariable::FieldVariable<FunctionSpaceType,3>>(this->functionSpace_->fieldVariable("general"));
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
