#include "data_management/mechanics.h"

namespace Data
{

template<typename FunctionSpaceType>
Mechanics<FunctionSpaceType>::
Mechanics(DihuContext context) : Data<FunctionSpaceType>(context)
{
}

template<typename FunctionSpaceType>
void Mechanics<FunctionSpaceType>::
initialize()
{
  Data<FunctionSpaceType>::initialize();
}

template<typename FunctionSpaceType>
void Mechanics<FunctionSpaceType>::
reset()
{
  // set initalize_ to false
  Data<FunctionSpaceType>::reset();

  // deallocate Petsc matrices
  this->stiffnessMatrix_ = nullptr;
}

template<typename FunctionSpaceType>
void Mechanics<FunctionSpaceType>::
setExternalForcesRightHandSide(std::shared_ptr<FieldVariableType> rhs)
{
  this->rhs_ = rhs;

  // multiply with -1
  PetscErrorCode ierr;
  ierr = VecScale(this->rhs_->valuesGlobal(), -1); CHKERRV(ierr);
}

template<typename FunctionSpaceType>
void Mechanics<FunctionSpaceType>::
createPetscObjects()
{
  assert(this->functionSpace_);

  // get the partitioning from the function space
  std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> meshPartition = this->functionSpace_->meshPartition();

  // create field variables on local partition
  this->displacements_ = this->functionSpace_->template createFieldVariable<FunctionSpaceType::dim()>("displacements");

  // create PETSc matrix object

  // PETSc MatCreateAIJ parameters
  int diagonalNonZeros = 3;   // number of nonzeros per row in DIAGONAL portion of local submatrix (same value is used for all local rows)
  int offdiagonalNonZeros = 0;   //  number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix (same value is used for all local rows)

  getPetscMemoryParameters(diagonalNonZeros, offdiagonalNonZeros);

  LOG(DEBUG) << "d=" << this->functionSpace_->dimension()
    << ", number of diagonal non-zeros: " << diagonalNonZeros << ", number of off-diagonal non-zeros: " <<offdiagonalNonZeros;

  const int nComponents = 1;
  const int D = FunctionSpaceType::dim();
  std::array<Mat,D*D> submatrices;

  for (int j = 0; j < D; j++)
  {
    for (int i = 0; i < D; i++)
    {
      std::stringstream name;
      name << "stiffnessMatrix_" << j << "," << i;

      stiffnessMatrixComponents_[j*D + i] = std::make_shared<PartitionedPetscMat<FunctionSpaceType>>(meshPartition, nComponents, diagonalNonZeros, offdiagonalNonZeros, name.str());
      submatrices[j*D + i] = stiffnessMatrixComponents_[j*D + i]->valuesGlobal();
    }
  }

  // create nested matrix
  PetscErrorCode ierr;
  ierr = MatCreateNest(meshPartition->mpiCommunicator(),
                       D, NULL, D, NULL, submatrices.data(), &this->stiffnessMatrix_); CHKERRV(ierr);
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,FunctionSpaceType::dim()>> Mechanics<FunctionSpaceType>::
displacements()
{
  return this->displacements_;
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,FunctionSpaceType::dim()>> Mechanics<FunctionSpaceType>::
rightHandSide()
{
  return this->rightHandSide_;
}

template<typename FunctionSpaceType>
Mat &Mechanics<FunctionSpaceType>::
stiffnessMatrix()
{
  return this->stiffnessMatrix_;
}

template<typename FunctionSpaceType>
void Mechanics<FunctionSpaceType>::
print()
{
  if (!VLOG_IS_ON(4))
    return;

  VLOG(4) << "======================";
  VLOG(4) << *this->displacements_;
  VLOG(4) << "======================";
}

template<typename FunctionSpaceType>
typename Mechanics<FunctionSpaceType>::OutputFieldVariables Mechanics<FunctionSpaceType>::
getOutputFieldVariables()
{
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> geometryField
    = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,3>>(this->functionSpace_->geometryField());

  return OutputFieldVariables(
    geometryField,
    displacements_,
    rhs_
  );
}

} // namespace Data

#include "data_management/mechanics.tpp"
