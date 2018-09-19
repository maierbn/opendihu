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
Mat &ModelOrderReduction<FullFunctionSpaceType>::
basis()
{    
  return this->basis_; 
}
  
template<typename FullFunctionSpaceType>  
Mat &ModelOrderReduction<FullFunctionSpaceType>::
basisTransp()
{    
  return this->basisTransp_; 
}
  
//template<typename FullFunctionSpaceType>  
//void ModelOrderReduction::setBasis()
//{
//}
  
template<typename FullFunctionSpaceType>  
Mat &ModelOrderReduction<FullFunctionSpaceType>::
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
  
template<typename FullFunctionSpaceType>
void ModelOrderReduction<FullFunctionSpaceType>::
createPetscObjects()
{
  LOG(TRACE) << "ModelOrderReduction::createPetscObjects()";
  
  // create field variables on local partition
  this->redSolution_ = this->functionSpace_->createGenericFieldVariable("redSolution");
  this->redIncrement_=this->functionSpace_->createGenericFieldVariable("redIncrement");
  
  // get the partitioning from the function space
  std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> meshPartition_row = this->fullFunctionSpace_->meshPartition();
  std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> meshPartition_col = this->functionSpace_->meshPartition();
  
  //this-> matrix = ... PetscPartitionedMat <FullFunctionSpaceType,FunctionSpace::Generic>
}
  
//template<typename FunctionSpcaeType>
//void initializeRedSysMatrix(Mat &A_R)
//{
  //to be implemented
//}

} //namespace
