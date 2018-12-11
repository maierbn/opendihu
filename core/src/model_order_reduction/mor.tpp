#include "model_order_reduction/mor.h"
#include "data_management/data.h"
//#include <petscmat.h>
#include <array>

namespace ModelOrderReduction
{

template<typename FullFunctionSpaceType>
MORBase<FullFunctionSpaceType>::
MORBase(DihuContext context):
//dataMOR_(std::make_shared<DataMOR>(context)),
initialized_(false)
{
  this->dataMOR_ = std::make_shared <DataMOR>(context);
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
  assert(dataMOR_);
  Mat &basisTransp=this->dataMOR_->basisTransp()->valuesGlobal();
  
  PetscErrorCode ierr;
  ierr=MatShift(basisTransp, 1); CHKERRV(ierr); //identitty matrix to check
  
  Mat &basis=this->dataMOR_->basis()->valuesGlobal();
  
  ierr=MatShift(basis, 1); CHKERRV(ierr); //identitty matrix to check
  
  /*
  PetscInt mat_sz_1, mat_sz_2;
  MatGetSize(basisTransp,&mat_sz_1,&mat_sz_2);
  LOG(DEBUG) << "mat_sz_1: " << mat_sz_1 << "mat_sz_2" << mat_sz_2 << "==============";
  */
  
  //to be implemented
  
}

template<typename FullFunctionSpaceType>
Data::ModelOrderReduction<FullFunctionSpaceType> &MORBase<FullFunctionSpaceType>::
dataMOR()
{
  return *dataMOR_;
}

template<typename FullFunctionSpaceType>
void MORBase<FullFunctionSpaceType>::
initialize()
{  
  if (initialized_)
    return;
  
  LOG(TRACE) << "MORBase::initialize()";
  
  dataMOR_->initialize();
  
  setBasis();
  
  initialized_=true;
  
}

template<typename FullFunctionSpaceType>
void MORBase<FullFunctionSpaceType>::
setRedSysMatrix(Mat &A, Mat &A_R)
{
  //to be implemented
  
  //MatGetRow  //(1) to get each row of the V_k^T or V_k.
  //VecGetValues //(2-1-a) To extract only the required part of vector. Output is probably an array not petsc vector. Indirect indexing is considered. It is local.
  //VecCreatMPIWithArray  //(2-2-a)creating petsc vec from array
  //VecMdot //(3) Vector dot products  cor computing V_k^T A. We do not need to stor the required parts of V_k and V_k^T as new matrices.
  
  //instead of (2-a) one can use the following:
  //ISCreatGeneral //(2-1-b)creat the index set for the indices that we require.
  //VecGetSubVector //(2-2-b) uses the above index set
  
}
template<typename FullFunctionSpaceType>
void MORBase<FullFunctionSpaceType>::
MatMultReduced(Mat mat,Vec x,Vec y)
{
  PetscErrorCode ierr;
  
  int vec_sz,mat_sz_1,mat_sz_2;
  VecGetSize(x,&vec_sz);
  MatGetSize(mat,&mat_sz_1,&mat_sz_2);
  
  if(mat_sz_2==vec_sz)
  {
    ierr=MatMult(mat,x,y);
    CHKERRV(ierr);
  }
  else
  {
    LOG(ERROR) << "MatMultReduced be done!";
  }
}

} //namespace
