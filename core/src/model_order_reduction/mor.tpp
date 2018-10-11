#include "model_order_reduction/mor.h"
#include "data_management/data.h"

#include <array>

namespace ModelOrderReduction
{

template<typename FullFunctionSpaceType>
MORBase<FullFunctionSpaceType>::
MORBase(DihuContext context):
  data_(std::make_shared<Data>(context)), initialized_(false)
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
  Mat &basisTransp=this->data_->basisTransp()->valuesGlobal();
  
  PetscErrorCode ierr;
  ierr=MatShift(basisTransp, 1); CHKERRV(ierr); //identitty matrix to check
  
  Mat &basis=this->data_->basis()->valuesGlobal();
  
  ierr=MatShift(basis, 1); CHKERRV(ierr); //identitty matrix to check
  
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
  
  setBasis();
  
  initialized_=true;
  
}

template<typename FullFunctionSpaceType>
void MORBase<FullFunctionSpaceType>::
setRedSysMatrix(Mat &A, Mat &A_R)
{
  //to be implemented
}

} //namespace
