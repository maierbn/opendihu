#pragma once

#include <petscmat.h>
#include "easylogging++.h"

#include "data_management/model_order_reduction.h"


namespace Data
{

template<typename FunctionSpaceType>  
ModelOrderReduction<FunctionSpaceType>::
ModelOrderReduction(DihuContext context): Data<FunctionSpaceType>(context)
{    
}

template<typename FunctionSpaceType>  
ModelOrderReduction<FunctionSpaceType>::
~ModelOrderReduction()
{    
}

template<typename FunctionSpaceType>  
Mat &ModelOrderReduction<FunctionSpaceType>::
basis()
{    
  return this->basis_; 
}
  
template<typename FunctionSpaceType>  
Mat &ModelOrderReduction<FunctionSpaceType>::
basisTransp()
{    
  return this->basisTransp_; 
}
  
//template<typename FunctionSpaceType>  
//void ModelOrderReduction::setBasis()
//{
//}
  
template<typename FunctionSpaceType>  
Mat &ModelOrderReduction<FunctionSpaceType>::
redSysMatrix()
{    
  return this->redSysMatrix_; 
} 
 
template<typename FunctionSpaceType>  
void ModelOrderReduction<FunctionSpaceType>::
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
  
template<typename FunctionSpaceType>
void ModelOrderReduction<FunctionSpaceType>::
createPetscObjects()
{
  LOG(TRACE) << "ModelOrderReduction::createPetscObjects()";
  
  // get the partitioning from the function space
  std::shared_ptr<Partition::MeshPartition<FunctionSpaceType>> meshPartition = this->functionSpace_->meshPartition();
    
  // create field variables on local partition
  this->redSolution_ = this->functionSpace_->template createFieldVariable<1>("redSolution"); 
}
  
//template<typename FunctionSpcaeType>
//void initializeRedSysMatrix(Mat &A_R)
//{
  //to be implemented
//}

} //namespace
