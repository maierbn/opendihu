#include "data_management/time_stepping/time_stepping_reduced.h"

#include <petscmat.h>
#include "easylogging++.h"

namespace Data
{

template<typename functionSpaceRowsType>  
TimeSteppingReduced<functionSpaceRowsType>::
TimeSteppingReduced(DihuContext context):
ModelOrderReduction<functionSpaceRowsType>(context),
TimeStepping<::FunctionSpace::Generic,1>(context)
{    
}

template<typename functionSpaceRowsType>  
TimeSteppingReduced<functionSpaceRowsType>::
~TimeSteppingReduced()
{    
}

template<typename functionSpaceRowsType> 
void TimeSteppingReduced<functionSpaceRowsType>::
setFunctionSpaceRows(std::shared_ptr<functionSpaceRowsType> functionSpace)
{
  this->functionSpaceRows_=functionSpace;
}

template<typename functionSpaceRowsType>  
std::shared_ptr<PartitionedPetscMat<functionSpaceRowsType,::FunctionSpace::Generic>> 
&TimeSteppingReduced<functionSpaceRowsType>::
basis()
{    
  return this->basis_; 
}
  
template<typename functionSpaceRowsType>  
std::shared_ptr<PartitionedPetscMat<::FunctionSpace::Generic,functionSpaceRowsType>> 
&TimeSteppingReduced<functionSpaceRowsType>::
basisTransp()
{    
  return this->basisTransp_; 
}
  
//template<typename functionSpaceRowsType>  
//void TimeSteppingReduced::setBasis()
//{
//}
  
template<typename functionSpaceRowsType>  
std::shared_ptr<PartitionedPetscMat<::FunctionSpace::Generic,::FunctionSpace::Generic>>
&TimeSteppingReduced<functionSpaceRowsType>::
redSysMatrix()
{    
  return this->redSysMatrix_; 
} 
 
template<typename functionSpaceRowsType>  
void TimeSteppingReduced<functionSpaceRowsType>::
initialize()
{
  if (!this->initialized_)
  {
    this->createPetscObjects();
    this->initialized_ = true;
  }
  else
  {
    LOG(WARNING) << "TimeSteppingReduced::Initialize(), TimeSteppingReduced is already assigned";
  }
}

//! the reduced solution
template<typename functionSpaceRowsType>  
std::shared_ptr<typename TimeSteppingReduced<functionSpaceRowsType>::FieldVariableType> &
TimeSteppingReduced<functionSpaceRowsType>::
redSolution()
{
  return redSolution_;
}
  
//! The reduced order increment
template<typename functionSpaceRowsType>  
std::shared_ptr<typename TimeSteppingReduced<functionSpaceRowsType>::FieldVariableType> &
TimeSteppingReduced<functionSpaceRowsType>::
redIncrement()
{
  return redIncrement_;
}

template<typename functionSpaceRowsType>
void TimeSteppingReduced<functionSpaceRowsType>::
createPetscObjects()
{
  LOG(TRACE) << "TimeSteppingReduced::createPetscObjects()";
  
  // create field variables on local partition
  const int nComponents = 1;
  this->redSolution_ = this->functionSpace_->template createFieldVariable<nComponents>("redSolution");
  this->redIncrement_ = std::static_pointer_cast<FieldVariableType>(this->functionSpace_->createFieldVariable("redIncrement", 1));
  
  // get the partitioning from the function space
  std::shared_ptr<Partition::MeshPartition<functionSpaceRowsType>> 
  meshPartitionRows = this->functionSpaceRows_->meshPartition();
  
  std::shared_ptr<Partition::MeshPartition<::FunctionSpace::Generic>>
  meshPartitionColumns = this->functionSpace_->meshPartition();
  
  this->basis_ = std::make_shared<PartitionedPetscMat<::FunctionSpace::Generic,functionSpaceRowsType>>(
    meshPartitionRows, meshPartitionColumns, 1, "basis");
  this->basisTransp_ = std::make_shared<PartitionedPetscMat<::FunctionSpace::Generic,functionSpaceRowsType>>(
    meshPartitionColumns, meshPartitionRows, 1, "basisTransp");
  this->redSysMatrix_=std::make_shared<PartitionedPetscMat<::FunctionSpace::Generic,::FunctionSpace::Generic>>(
    meshPartitionColumns, meshPartitionColumns, 1, "redSysMatrix");
}

} //namespace
