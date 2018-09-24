#include "data_management/model_order_reduction.h"

#include <petscmat.h>
#include "easylogging++.h"

namespace Data
{

template<typename FullFunctionSpaceType>  
ModelOrderReduction<FullFunctionSpaceType>::
ModelOrderReduction(DihuContext context): Data<FullFunctionSpaceType>(context)
{    
}

template<typename FullFunctionSpaceType>  
ModelOrderReduction<FullFunctionSpaceType>::
~ModelOrderReduction()
{    
}

template<typename FullFunctionSpaceType> 
void ModelOrderReduction<FullFunctionSpaceType>::
setFullFunctionSpace(std::shared_ptr<FullFunctionSpaceType> functionSpace)
{
  this->fullFunctionSpace_=functionSpace;
}

template<typename FullFunctionSpaceType>  
std::shared_ptr<PartitionedPetscMat<FullFunctionSpaceType,::FunctionSpace::Generic>> 
&ModelOrderReduction<FullFunctionSpaceType>::
basis()
{    
  return this->basis_; 
}
  
template<typename FullFunctionSpaceType>  
std::shared_ptr<PartitionedPetscMat<::FunctionSpace::Generic,FullFunctionSpaceType>> 
&ModelOrderReduction<FullFunctionSpaceType>::
basisTransp()
{    
  return this->basisTransp_; 
}
  
//template<typename FullFunctionSpaceType>  
//void ModelOrderReduction::setBasis()
//{
//}
  
template<typename FullFunctionSpaceType>  
std::shared_ptr<PartitionedPetscMat<::FunctionSpace::Generic,::FunctionSpace::Generic>>
&ModelOrderReduction<FullFunctionSpaceType>::
redSysMatrix()
{    
  return this->redSysMatrix_; 
} 
 
template<typename FullFunctionSpaceType>  
void ModelOrderReduction<FullFunctionSpaceType>::
initialize()
{
  if (!this->initialized_)
  {
    this->createPetscObjects();
    this->initialized_ = true;
  }
  else
  {
    LOG(WARNING) << "ModelOrderReduction::Initialize(), ModelOrderReduction is already assigned";
  }
}

//! the reduced solution
template<typename FullFunctionSpaceType>  
std::shared_ptr<typename ModelOrderReduction<FullFunctionSpaceType>::FieldVariableType> &
ModelOrderReduction<FullFunctionSpaceType>::
redSolution()
{
  return redSolution_;
}
  
//! The reduced order increment
template<typename FullFunctionSpaceType>  
std::shared_ptr<typename ModelOrderReduction<FullFunctionSpaceType>::FieldVariableType> &
ModelOrderReduction<FullFunctionSpaceType>::
redIncrement()
{
  return redIncrement_;
}

template<typename FullFunctionSpaceType>
void ModelOrderReduction<FullFunctionSpaceType>::
createPetscObjects()
{
  LOG(TRACE) << "ModelOrderReduction::createPetscObjects()";
  
  // create field variables on local partition
  const int nComponents = 1;
  this->redSolution_ = this->functionSpace_->template createFieldVariable<nComponents>("redSolution");
  this->redIncrement_ = std::static_pointer_cast<FieldVariableType>(this->functionSpace_->createFieldVariable("redIncrement", 1));
  
  // get the partitioning from the function space
  std::shared_ptr<Partition::MeshPartition<FullFunctionSpaceType>> 
  meshPartitionRows = this->fullFunctionSpace_->meshPartition();
  
  std::shared_ptr<Partition::MeshPartition<::FunctionSpace::Generic>>
  meshPartitionColumns = this->functionSpace_->meshPartition();
  
  this->basis_ = std::make_shared<PartitionedPetscMat<::FunctionSpace::Generic,FullFunctionSpaceType>>(
    meshPartitionRows, meshPartitionColumns, 1, "basis");
  this->basisTransp_ = std::make_shared<PartitionedPetscMat<::FunctionSpace::Generic,FullFunctionSpaceType>>(
    meshPartitionColumns, meshPartitionRows, 1, "basisTransp");
  this->redSysMatrix_=std::make_shared<PartitionedPetscMat<::FunctionSpace::Generic,::FunctionSpace::Generic>>(
    meshPartitionColumns, meshPartitionColumns, 1, "redSysMatrix");
}

} //namespace
