#include "data_management/time_stepping/time_stepping_implicit.h"

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

template<typename FunctionSpaceType,int nComponents>
void TimeSteppingImplicit<FunctionSpaceType,nComponents>::
createPetscObjects()
{
  LOG(DEBUG) << "TimeSteppingImplicit::createPetscObjects(" <<nComponents << ")";
  TimeStepping<FunctionSpaceType,nComponents>::createPetscObjects();

  this->boundaryConditionsRightHandSideSummand_ = this->functionSpace_->template createFieldVariable<nComponents>("boundaryConditionsRightHandSideSummand");
  this->systemRightHandSide_ = this->functionSpace_->template createFieldVariable<nComponents>("systemRightHandSide");

  this->debuggingName_ = "Implicit";
  VLOG(1) << "initial values " << *this->boundaryConditionsRightHandSideSummand_;
}

template<typename FunctionSpaceType,int nComponents>
std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>>
TimeSteppingImplicit<FunctionSpaceType,nComponents>::
systemMatrix()
{
  return this->systemMatrix_;
}

template<typename FunctionSpaceType,int nComponents>
std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>>
TimeSteppingImplicit<FunctionSpaceType,nComponents>::
integrationMatrixRightHandSide()
{
  return this->integrationMatrixRightHandSide_;
}

template<typename FunctionSpaceType,int nComponents>
void TimeSteppingImplicit<FunctionSpaceType,nComponents>::
initializeSystemMatrix(Mat &systemMatrix)
{
  // if the systemMatrix_ is already initialized do not initialize again
  if (this->systemMatrix_)
  {
    LOG(WARNING) << "Initialize system matrix again. Previous system matrix: " << *this->systemMatrix_;
  }

  // the PETSc matrix object is created outside by MatMatMult
  std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> partition = this->functionSpace_->meshPartition();
  this->systemMatrix_ = std::make_shared<PartitionedPetscMat<FunctionSpaceType>>(partition, systemMatrix, "systemMatrix");
}

template<typename FunctionSpaceType,int nComponents>
void TimeSteppingImplicit<FunctionSpaceType,nComponents>::
initializeIntegrationMatrixRightHandSide(Mat &integrationMatrix)
{
  // if the integrationMatrix_ is already initialized do not initialize again
  if (this->integrationMatrixRightHandSide_)
    return;

  // the PETSc matrix object is created outside by MatMatMult
  std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> partition = this->functionSpace_->meshPartition();
  this->integrationMatrixRightHandSide_ = std::make_shared<PartitionedPetscMat<FunctionSpaceType>>(partition, integrationMatrix, "integrationMatrixRightHandSide");
}

template<typename FunctionSpaceType,int nComponents>
void TimeSteppingImplicit<FunctionSpaceType,nComponents>::
initializeMatrix(Mat &matrixIn, std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> matrixOut, std::string name)
{
  // if the matrix is already initialized do not initialize again
  if (matrixOut)
    return;

  // the PETSc matrix object is created somewhere else
  std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> partition = this->functionSpace_->meshPartition();
  matrixOut = std::make_shared<PartitionedPetscMat<FunctionSpaceType>>(partition, matrixIn, name);
}
template<typename FunctionSpaceType,int nComponents>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> TimeSteppingImplicit<FunctionSpaceType,nComponents>::
boundaryConditionsRightHandSideSummand()
{
  return this->boundaryConditionsRightHandSideSummand_;
}

template<typename FunctionSpaceType,int nComponents>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> TimeSteppingImplicit<FunctionSpaceType,nComponents>::
systemRightHandSide()
{
  return this->systemRightHandSide_;
}

template<typename FunctionSpaceType,int nComponents>
void TimeSteppingImplicit<FunctionSpaceType,nComponents>::
print() // use override in stead of extending the parents' print output.This way "solution" is still in the end.
{
  if (!VLOG_IS_ON(4))
    return;

  TimeStepping<FunctionSpaceType,nComponents>::print();

  if (this->systemMatrix_)
    VLOG(4) << *this->systemMatrix_;
  VLOG(4) << *this->boundaryConditionsRightHandSideSummand_;
  VLOG(4) << "=======================================";
}

} // namespace
