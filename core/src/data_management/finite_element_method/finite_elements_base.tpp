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

template<typename FunctionSpaceType, int nComponents>
FiniteElementsBase<FunctionSpaceType,nComponents>::
FiniteElementsBase(DihuContext context) : Data<FunctionSpaceType>(context)
{
}

template<typename FunctionSpaceType, int nComponents>
FiniteElementsBase<FunctionSpaceType,nComponents>::
~FiniteElementsBase()
{
}

template<typename FunctionSpaceType, int nComponents>
void FiniteElementsBase<FunctionSpaceType,nComponents>::
initialize()
{
  LOG(DEBUG) << "Data::FiniteElementsBase::initialize";
  Data<FunctionSpaceType>::initialize();

  slotConnectorData_ = std::make_shared<SlotConnectorDataType>();
  slotConnectorData_->addFieldVariable(this->solution());

  // parse slot names of the additional field variables
  std::string slotName = this->context_.getPythonConfig().getOptionString("slotName", "");
  slotConnectorData_->slotNames.push_back(slotName);
}
template<typename FunctionSpaceType, int nComponents>
void FiniteElementsBase<FunctionSpaceType,nComponents>::
reset()
{
  LOG(DEBUG) << "Data::FiniteElementsBase::reset";
  // set initalize_ to false
  Data<FunctionSpaceType>::reset();

  // deallocate Petsc matrices
  this->stiffnessMatrixWithoutBc_ = nullptr;
  this->stiffnessMatrix_ = nullptr;
  LOG(DEBUG) << "stiffnessMatrix_ set to nullptr";
}

template<typename FunctionSpaceType, int nComponents>
void FiniteElementsBase<FunctionSpaceType,nComponents>::
getPetscMemoryParameters(int &nNonZerosDiagonal, int &nNonZerosOffdiagonal)
{
  const int D = FunctionSpaceType::dim();
  const int nDofsPerNode = FunctionSpace::FunctionSpaceBaseDim<1,typename FunctionSpaceType::BasisFunction>::nDofsPerNode();
  const int nDofsPerBasis = FunctionSpace::FunctionSpaceBaseDim<1,typename FunctionSpaceType::BasisFunction>::nDofsPerElement();
  const int nOverlaps = (nDofsPerBasis*2 - 1) * nDofsPerNode;   // number of nodes of 2 neighbouring 1D elements (=number of ansatz functions in support of center ansatz function)

  // due to PETSc storage nNonZerosDiagonal and nNonZerosOffdiagonal should be both set to the maximum number of non-zero entries per row

  switch (D)
  {
  case 1:
    nNonZerosDiagonal = nOverlaps;
    nNonZerosOffdiagonal = nNonZerosDiagonal;
    break;
  case 2:
    nNonZerosDiagonal = pow(nOverlaps, 2) + 16;   // because of boundary conditions there can be more entries, which are all zero, but stored as non-zero
    nNonZerosOffdiagonal = nNonZerosDiagonal;
    break;
  case 3:
    nNonZerosDiagonal = pow(nOverlaps, 3);
    nNonZerosOffdiagonal = nNonZerosDiagonal;
    break;
  };
}

template<typename FunctionSpaceType, int nComponents>
void FiniteElementsBase<FunctionSpaceType,nComponents>::
createPetscObjects()
{
  LOG(TRACE) << "FiniteElements::createPetscObjects";

  assert(this->functionSpace_);

  // get the partitioning from the function space
  std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> meshPartition = this->functionSpace_->meshPartition();
  
  // create field variables on local partition
  this->rhs_ = this->functionSpace_->template createFieldVariable<nComponents>("rightHandSide");
  this->solution_ = this->functionSpace_->template createFieldVariable<nComponents>("solution");
  this->negativeRhsNeumannBoundaryConditions_ = this->functionSpace_->template createFieldVariable<nComponents>("zero");

  // create PETSc matrix object

  // PETSc MatCreateAIJ parameters
  int nNonZerosDiagonal = 3;   // number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
  int nNonZerosOffdiagonal = 0;   //  number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows)

  getPetscMemoryParameters(nNonZerosDiagonal, nNonZerosOffdiagonal);

  LOG(DEBUG) << "d=" << this->functionSpace_->dimension()
    << ", number of diagonal non-zeros: " << nNonZerosDiagonal << ", number of off-diagonal non-zeros: " <<nNonZerosOffdiagonal;

  LOG(DEBUG) << "create new stiffnessMatrix";
  this->stiffnessMatrix_ = std::make_shared<PartitionedPetscMat<FunctionSpaceType>>(meshPartition, nComponents, nNonZerosDiagonal, nNonZerosOffdiagonal, "stiffnessMatrix");
  this->stiffnessMatrixWithoutBc_ = std::make_shared<PartitionedPetscMat<FunctionSpaceType>>(meshPartition, nComponents, nNonZerosDiagonal, nNonZerosOffdiagonal, "stiffnessMatrixWithoutBc");
}

template<typename FunctionSpaceType, int nComponents>
std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> FiniteElementsBase<FunctionSpaceType,nComponents>::
stiffnessMatrix()
{
  return this->stiffnessMatrix_;
}

template<typename FunctionSpaceType, int nComponents>
std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> FiniteElementsBase<FunctionSpaceType,nComponents>::
stiffnessMatrixWithoutBc()
{
  return this->stiffnessMatrixWithoutBc_;
}

template<typename FunctionSpaceType, int nComponents>
std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> FiniteElementsBase<FunctionSpaceType,nComponents>::
massMatrix()
{
  return this->massMatrix_;
}

template<typename FunctionSpaceType, int nComponents>
std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> FiniteElementsBase<FunctionSpaceType,nComponents>::
inverseLumpedMassMatrix()
{
  return this->inverseLumpedMassMatrix_;
}

template<typename FunctionSpaceType, int nComponents>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> FiniteElementsBase<FunctionSpaceType,nComponents>::
rightHandSide()
{
  return this->rhs_;
}

template<typename FunctionSpaceType, int nComponents>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> FiniteElementsBase<FunctionSpaceType,nComponents>::
negativeRightHandSideNeumannBoundaryConditions()
{
  return this->negativeRhsNeumannBoundaryConditions_;
}

template<typename FunctionSpaceType, int nComponents>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> FiniteElementsBase<FunctionSpaceType,nComponents>::
solution()
{
  return this->solution_;
}

//! set the solution variable if it is initialized externally, such as in a timestepping scheme
template<typename FunctionSpaceType, int nComponents>
void FiniteElementsBase<FunctionSpaceType,nComponents>::
setSolutionVariable(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> solution)
{
  // this will be called by the time stepping scheme after initialize()
  this->solution_ = solution;
}

template<typename FunctionSpaceType, int nComponents>
void FiniteElementsBase<FunctionSpaceType,nComponents>::
setNegativeRightHandSideNeumannBoundaryConditions(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> negativeRightHandSideNeumannBoundaryConditions)
{
  this->negativeRhsNeumannBoundaryConditions_ = negativeRightHandSideNeumannBoundaryConditions;
}

template<typename FunctionSpaceType, int nComponents>
std::shared_ptr<typename FiniteElementsBase<FunctionSpaceType,nComponents>::SlotConnectorDataType>
FiniteElementsBase<FunctionSpaceType,nComponents>::
getSlotConnectorData()
{
  return this->slotConnectorData_;
}

template<typename FunctionSpaceType, int nComponents>
void FiniteElementsBase<FunctionSpaceType,nComponents>::
print()
{
  if (!VLOG_IS_ON(4))
    return;

  VLOG(4) << "======================";
  VLOG(4) << "nComponents: " << nComponents;
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

template<typename FunctionSpaceType, int nComponents>
void FiniteElementsBase<FunctionSpaceType,nComponents>::
initializeMassMatrix()
{
  // if the massMatrix is already initialized do not initialize again
  if (this->massMatrix_)
    return;

  // create PETSc matrix object

  // PETSc MatCreateAIJ parameters
  int nNonZerosDiagonal = 3;   // number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
  int nNonZerosOffdiagonal = 0;   //  number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows)

  getPetscMemoryParameters(nNonZerosDiagonal, nNonZerosOffdiagonal);

  std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> partition = this->functionSpace_->meshPartition();
  this->massMatrix_ = std::make_shared<PartitionedPetscMat<FunctionSpaceType>>(partition, nComponents, nNonZerosDiagonal, nNonZerosOffdiagonal, "massMatrix");
}

template<typename FunctionSpaceType, int nComponents>
void FiniteElementsBase<FunctionSpaceType,nComponents>::
initializeInverseLumpedMassMatrix()
{
  // if the inverseLumpedMassMatrix_ is already initialized do not initialize again
  if (this->inverseLumpedMassMatrix_)
    return;

  // create PETSc matrix object

  // PETSc MatCreateAIJ parameters
  int nNonZerosDiagonal = 3;   // number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
  int nNonZerosOffdiagonal = 0;   //  number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows)

  getPetscMemoryParameters(nNonZerosDiagonal, nNonZerosOffdiagonal);

  assert(this->functionSpace_);
  std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> partition = this->functionSpace_->meshPartition();
  this->inverseLumpedMassMatrix_ = std::make_shared<PartitionedPetscMat<FunctionSpaceType>>(partition, nComponents, nNonZerosDiagonal, nNonZerosOffdiagonal, "inverseLumpedMassMatrix");
}

template<typename FunctionSpaceType, int nComponents>
typename FiniteElementsBase<FunctionSpaceType,nComponents>::FieldVariablesForOutputWriter FiniteElementsBase<FunctionSpaceType,nComponents>::
getFieldVariablesForOutputWriter()
{
  // these field variables will be written to output files
  assert(this->functionSpace_);
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> geometryField
    = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,3>>(this->functionSpace_->geometryField());

  return FieldVariablesForOutputWriter(
    geometryField,
    solution_,
    rhs_,
    negativeRhsNeumannBoundaryConditions_
  );
}


} // namespace Data
