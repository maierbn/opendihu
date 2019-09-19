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
redSystemMatrix()
{    
  return this->redSystemMatrix_; 
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

template<typename FunctionSpaceRows>
void ModelOrderReduction<FunctionSpaceRows>::
createPetscObjects()
{
  LOG(TRACE) << "ModelOrderReduction::createPetscObjects()";
  
  assert(this->functionSpace_);
  assert(this->functionSpaceRows_);

  // get the partitioning from the function space
  std::shared_ptr<Partition::MeshPartition<FunctionSpaceRows>> 
  meshPartitionRows = this->functionSpaceRows_->meshPartition();
  
  std::shared_ptr<Partition::MeshPartition<::FunctionSpace::Generic>>
  meshPartitionColumns = this->functionSpace_->meshPartition();
  
  this->basis_ = std::make_shared<PartitionedPetscMat<FunctionSpaceRows,::FunctionSpace::Generic>>(
    meshPartitionRows, meshPartitionColumns, 1, "basis");
  this->basisTransp_ = std::make_shared<PartitionedPetscMat<::FunctionSpace::Generic,FunctionSpaceRows>>(
    meshPartitionColumns, meshPartitionRows, 1, "basisTransp");
  this->redSystemMatrix_=std::make_shared<PartitionedPetscMat<::FunctionSpace::Generic,::FunctionSpace::Generic>>(
    meshPartitionColumns, meshPartitionColumns, 1, "redSysMatrix");
}

} //namespace
