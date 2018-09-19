#include "model_order_reduction/mor.h"
#include "data_management/data.h"

#include <array>

namespace ModelOrderReduction
{

template<typename FullFunctionSpaceType>
MORBase<FullFunctionSpaceType>::
MORBase(DihuContext context):
initialized_(false)
{ 
}

template<typename FullFunctionSpaceType>
MORBase<FullFunctionSpaceType>
::~MORBase()
{    
}

template<typename FullFunctionSpaceType>
void MORBase<FullFunctionSpaceType>::
setBasis()
{
    //to be implemented
}

template<typename FullFunctionSpaceType>
Data::ModelOrderReduction<FullFunctionSpaceType> &MORBase<FullFunctionSpaceType>::
data()
{
  return *data_;
}

template<typename FullFunctionSpaceType>
void MORBase<FullFunctionSpaceType>::
initialize()
{  
  if (initialized_)
    return;
  
  LOG(TRACE) << "MORBase::initialize()";
  
}

template<typename FullFunctionSpaceType>
void MORBase<FullFunctionSpaceType>::
setRedSysMatrix(Mat &A, Mat &A_R)
{
  //to be implemented
}

} //namespace
