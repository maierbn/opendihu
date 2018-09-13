#include "model_order_reduction/mor.h"

namespace ModelOrderReduction
{

template<typename FunctionSpaceType>
MORBase<FunctionSpaceType>::
MORBase(DihuContext context):
initialized_(false)
{    
}

template<typename FunctionSpaceType>
MORBase<FunctionSpaceType>
::~MORBase()
{    
}

template<typename FunctionSpaceType>
void MORBase<FunctionSpaceType>::
setBasis()
{
    //to be implemented
}

template<typename FunctionSpaceType>
Data::ModelOrderReduction<FunctionSpaceType> &MORBase<FunctionSpaceType>::
data()
{
  return *data_;
}

template<typename FunctionSpaceType>
void MORBase<FunctionSpaceType>::
initialize()
{   
}

template<typename FunctionSpaceType>
void MORBase<FunctionSpaceType>::
setRedSysMatrix(Mat &A, Mat &A_R)
{
  //to be implemented
}

} //namespace
