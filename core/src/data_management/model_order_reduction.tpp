#include "data_management/model_order_reduction.h"

#include <petscmat.h>
#include "easylogging++.h"

namespace Data
{

template<typename FunctionSpaceRows>  
ModelOrderReduction<FunctionSpaceRows>::
ModelOrderReduction(DihuContext context): Data<FunctionSpaceRows>(context)
{    
}

template<typename FunctionSpaceRows>  
ModelOrderReduction<FunctionSpaceRows>::
~ModelOrderReduction()
{    
}

template<typename FunctionSpaceRows> 
void ModelOrderReduction<FunctionSpaceRows>::
setFunctionSpaceRows(std::shared_ptr<FunctionSpaceRows> functionSpace)
{
  this->functionSpaceRows_=functionSpace;
}

template<typename FunctionSpaceRows>  
std::shared_ptr<PartitionedPetscMat<FunctionSpaceRows,::FunctionSpace::Generic>> 
&ModelOrderReduction<FunctionSpaceRows>::
basis()
{    
  return this->basis_; 
}
  
template<typename FunctionSpaceRows>  
std::shared_ptr<PartitionedPetscMat<::FunctionSpace::Generic,FunctionSpaceRows>> 
&ModelOrderReduction<FunctionSpaceRows>::
basisTransp()
{    
  return this->basisTransp_; 
}
  
//template<typename FunctionSpaceRows>  
//void ModelOrderReduction::setBasis()
//{
//}
  
template<typename FunctionSpaceRows>  
std::shared_ptr<PartitionedPetscMat<::FunctionSpace::Generic,::FunctionSpace::Generic>>
&ModelOrderReduction<FunctionSpaceRows>::
redSysMatrix()
{    
  return this->redSysMatrix_; 
} 
 
template<typename FunctionSpaceRows>  
void ModelOrderReduction<FunctionSpaceRows>::
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
template<typename FunctionSpaceRows>  
std::shared_ptr<typename ModelOrderReduction<FunctionSpaceRows>::FieldVariableType> &
ModelOrderReduction<FunctionSpaceRows>::
redSolution()
{
  return redSolution_;
}
  
//! The reduced order increment
template<typename FunctionSpaceRows>  
std::shared_ptr<typename ModelOrderReduction<FunctionSpaceRows>::FieldVariableType> &
ModelOrderReduction<FunctionSpaceRows>::
redIncrement()
{
  return redIncrement_;
}

template<typename FunctionSpaceRows>
void ModelOrderReduction<FunctionSpaceRows>::
createPetscObjects()
{
  LOG(TRACE) << "ModelOrderReduction::createPetscObjects()";
  
  // create field variables on local partition
  const int nComponents = 1;
  this->redSolution_ = this->functionSpace_->template createFieldVariable<nComponents>("redSolution");
  this->redIncrement_ = std::static_pointer_cast<FieldVariableType>(this->functionSpace_->createFieldVariable("redIncrement", 1));
  
  // get the partitioning from the function space
  std::shared_ptr<Partition::MeshPartition<FunctionSpaceRows>> 
  meshPartitionRows = this->functionSpaceRows_->meshPartition();
  
  std::shared_ptr<Partition::MeshPartition<::FunctionSpace::Generic>>
  meshPartitionColumns = this->functionSpace_->meshPartition();
  
  this->basis_ = std::make_shared<PartitionedPetscMat<::FunctionSpace::Generic,FunctionSpaceRows>>(
    meshPartitionRows, meshPartitionColumns, 1, "basis");
  this->basisTransp_ = std::make_shared<PartitionedPetscMat<::FunctionSpace::Generic,FunctionSpaceRows>>(
    meshPartitionColumns, meshPartitionRows, 1, "basisTransp");
  this->redSysMatrix_=std::make_shared<PartitionedPetscMat<::FunctionSpace::Generic,::FunctionSpace::Generic>>(
    meshPartitionColumns, meshPartitionColumns, 1, "redSysMatrix");
}

} //namespace
