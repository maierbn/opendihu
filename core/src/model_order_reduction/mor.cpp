#include "model_order_reduction/mor.h"

namespace ModelOrderReduction
{

  MORBase::MORBase(DihuContext context, int n, int k){
    
  }
  
  Mat &MORBase::basis(){
    
    return this->basis_; 
  }
  
  Mat &MORBase::basisTransp(){
    
    return this->basisTransp_; 
  }
  
  Mat &MORBase::redSysMatrix(){
    
    return this->redSysMatrix_; 
  }
  
  void MORBase::setBasis(){
    //to be implemented
  }
  
  void MORBase::createPetscObjects(){
    //to be implemented
  }
  
  void MORBase::setRedSysMatrix(Mat A, Mat A_R){
    //to be implemented
  }

} //namespace
