#include "model_order_reduction/mor.h"
#include "data_management/data.h"
//#include <petscmat.h>
#include <array>

namespace ModelOrderReduction
{

template<typename FunctionSpaceRowsType>
MORBase<FunctionSpaceRowsType>::
MORBase(DihuContext context):
//dataMOR_(std::make_shared<DataMOR>(context)),
initialized_(false)
{
  this->dataMOR_ = std::make_shared <DataMOR>(context);
}

template<typename FunctionSpaceRowsType>
MORBase<FunctionSpaceRowsType>
::~MORBase()
{    
}

template<typename FunctionSpaceRowsType>
void MORBase<FunctionSpaceRowsType>::
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

template<typename FunctionSpaceRowsType>
Data::ModelOrderReduction<FunctionSpaceRowsType> &MORBase<FunctionSpaceRowsType>::
dataMOR()
{
  return *dataMOR_;
}

template<typename FunctionSpaceRowsType>
void MORBase<FunctionSpaceRowsType>::
initialize()
{  
  if (initialized_)
    return;
  
  LOG(TRACE) << "MORBase::initialize()";
  
  dataMOR_->initialize();
  
  setBasis();
  
  initialized_=true;
  
}

template<typename FunctionSpaceRowsType>
void MORBase<FunctionSpaceRowsType>::
setRedSysMatrix(Mat &A, Mat &A_R)
{
  //to be implemented
  
  //MatGetRow  //(1) to get each row of the V_k^T or V_k. Provides an array out of the row.
  //VecGetValues //(2-1-a) To extract only the required part of vector. Output is probably an array not petsc vector. Indirect indexing is considered. It is local.
  //VecCreateMPIWithArray  //(2-2-a)creating petsc vec from array
  //VecMdot //(3) Vector dot products  for computing V_k^T A. We do not need to stor the required parts of V_k and V_k^T as new matrices.
  
  //instead of (2-a) one can use the following:
  //ISCreatGeneral //(2-1-b)creat the index set for the indices that we require.
  //VecGetSubVector //(2-2-b) uses the above index set
  
}
template<typename FunctionSpaceRowsType>
void MORBase<FunctionSpaceRowsType>::
MatMultReduced(Mat mat,Vec x,Vec y)
{
  PetscErrorCode ierr;
  
  PetscInt vec_sz,mat_sz_1,mat_sz_2;
  VecGetSize(x,&vec_sz);
  MatGetSize(mat,&mat_sz_1,&mat_sz_2);
  
  LOG(TRACE) << "mat_sz_2: " << mat_sz_2 << " vec_sz: " << vec_sz;
  
  if(mat_sz_2==vec_sz)
  {
    ierr=MatMult(mat,x,y);
    CHKERRV(ierr);
  }
  else
  {
    LOG(TRACE) << "MatMultReduced";
    
    const PetscInt *cols;    
    const PetscScalar *mat_row;
    Vec mat_row_vec; //to allocate
    
    PetscScalar val[mat_sz_1];
    // to allocate ??
        
    for(int i=0; i<mat_sz_1; i++)
    {
      ierr=MatGetRow(mat,i,NULL,&cols,&mat_row);  CHKERRV(ierr);
      // use only the required part of the array      
      ierr=VecCreateSeqWithArray(PETSC_COMM_SELF,vec_sz,vec_sz,mat_row,&mat_row_vec); CHKERRV(ierr);
      ierr=VecTDot(x,mat_row_vec,&val[i]);  CHKERRV(ierr);
    }
    
    ierr=VecCreateSeqWithArray(PETSC_COMM_SELF,mat_sz_1,mat_sz_1,val,&y);
    //VecCreateMPIWithArray ???
    
    //LOG(ERROR) << "MORBase::MatMultReduced to be implemented!";
  }
}

template<typename FunctionSpaceRowsType>
void MORBase<FunctionSpaceRowsType>::
MatMultFull(Mat mat,Vec x,Vec y)
{
  PetscErrorCode ierr;
  
  PetscInt vec_sz,mat_sz_1,mat_sz_2;
  VecGetSize(y,&vec_sz);
  MatGetSize(mat,&mat_sz_1,&mat_sz_2);
  
  LOG(TRACE) << "mat_sz_2: " << mat_sz_2 << " vec_sz: " << vec_sz;
  
  if(mat_sz_1==vec_sz)
  {
    ierr=MatMult(mat,x,y);
    CHKERRV(ierr);
  }
  else
  {
    const PetscInt *cols;    
    const PetscScalar *mat_row;
    Vec mat_row_vec; //to allocate
    
    PetscScalar val[vec_sz];
    // to allocate ??
    
    for(int i=0; i<vec_sz; i++)
    {
      ierr=MatGetRow(mat,i,NULL,&cols,&mat_row);  CHKERRV(ierr);      
      ierr=VecCreateSeqWithArray(PETSC_COMM_SELF,mat_sz_2,mat_sz_2,mat_row,&mat_row_vec); CHKERRV(ierr);
      ierr=VecTDot(x,mat_row_vec,&val[i]);  CHKERRV(ierr);
    }
    
    ierr=VecCreateSeqWithArray(PETSC_COMM_SELF,vec_sz,vec_sz,val,&y);
    //VecCreateMPIWithArray ???
    //LOG(ERROR) << "MORBase::MatMultFull to be implemented!";
  }
}

} //namespace
