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
  
  
  //MatTransposeGetMat(Mat A,Mat *M)
  
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
/*
template<typename FunctionSpaceRowsType>
void MORBase<FunctionSpaceRowsType>::
setRedSystemMatrix(std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> ptr_systemMatrix, std::shared_ptr<PartitionedPetscMat<FunctionSpaceType>> ptr_redSystemMatrix)
{
  //to be implemented
  Mat &basis = this->dataMOR_->basis()->valuesGlobal();
  Mat &basisTransp = this->dataMOR_->basisTransp()->valuesGlobal(); 
  
  Mat systemMatrix=ptr_systemMatrix->valuesGlobal();
  Mat redSystemMatrix=ptr_redSystemMatrix->valuesGlobal();
  
  PetscErrorCode ierr; 
  
  PetscInt mat_sz_1, mat_sz_2;
  
  MatGetSize(basisTransp,&mat_sz_1,&mat_sz_2);
  PetscInt nrows_bs=mat_sz_1;
  PetscInt ncols_bs=mat_sz_2;
  
  MatGetSize(systemMatrix,&mat_sz_1,&mat_sz_2);
  PetscInt nrows_sys=mat_sz_1; //square matrix
  //PetscInt ncols_sys=mat_sz_2; 
  
  if(ncols_bs == nrows_sys)
  {
    //! Reduction of the system matrix in case of compatible row spaces of the system matrix and reduced basis    
    //D=A*B*C
    MatMatMatMult(basisTransp,systemMatrix,basis,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&redSystemMatrix); CHKERRV(ierr);      
  }
  else
  {
    
    ////const std::shared_ptr<FunctionSpaceType> functionSpaceReduced=this->dataMOR->FunctionSpace();  
    //std::shared_ptr<Partition::MeshPartition<::FunctionSpace::Generic>>
    //meshPartitionRows_=this->dataMOR_->basisTransp()->meshPartitionRows();
    
    std::shared_ptr<Partition::MeshPartition<ptr_systemMatrix.RowsFunctionSpace>>
    meshPartitionRows_=redSystemMatrix.meshPartitionRows();
    
    std::shared_ptr<Partition::MeshPartition<ptr_systemMatrix.RowsFunctionSpace>> 
    meshPartitionColumns_=ptr_systemMatrix.meshPartitionRows(); //compatible for multiplication  
    
    std::shared_ptr<PartitionedPetscMat<ptr_systemMatrix.RowsFunctionSpace,redSystemMatrix.RowsFunctionSpace>> basis_submatrix;
    std::shared_ptr<PartitionedPetscMat<::FunctionSpace::Generic,ptr_systemMatrix.RowsFunctionSpace>> basisTransp_submatrix;
    
    Mat &basis_sbm=basis_submatrix->valuesGlobal();
    Mat &basisTransp_sbm=basisTransp_submatrix->valuesGlobal();
    
    const PetscScalar vals_total[ncols_bs];
    const PetscScalar vals_sbm[nrows_sys];
    
    for( int row=0; row<nrows_bs;row++)
    {
      // to get each row of the V_k^T
      MatGetRow(basisTransp,row,Null,Null,&vals_total); CHKERRV(ierr);    
      MatSetValuesRow(basisTransp_sbm,row,vals_total); CHKERRV(ierr); //inconsistent sizes may work. Otherwise build vals_sbm.     
    }
    
    MatTransposeGetMat(basisTransp_sbm,&basis_sbm); CHKERRV(ierr);
    
    //D=A*B*C
    MatMatMatMult(basisTransp_sbm,systemMatrix,basis_sbm,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&redSystemMatrix); CHKERRV(ierr);
  }
  
}
*/
template<typename FunctionSpaceRowsType>
void MORBase<FunctionSpaceRowsType>::
MatMultReduced(Mat mat,Vec x,Vec y)
{
  PetscErrorCode ierr;
  //PetscInt vec_sz_tmp;
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
    
    PetscInt *idx;
    PetscMalloc1(vec_sz,&idx);
    for(PetscInt i=0; i<vec_sz; i++)
          idx[i]=i;
    
    PetscInt *idx_2;
    PetscMalloc1(mat_sz_1,&idx_2);
    for(PetscInt i=0; i<mat_sz_1; i++)
      idx_2[i]=i;
           
    const PetscScalar *mat_row;
    PetscScalar mat_row_sbv[vec_sz];
    
    Vec mat_row_vec;
    ierr=VecDuplicate(x,&mat_row_vec); CHKERRV(ierr);
    
    PetscScalar val[mat_sz_1];
       
    for(int i=0; i<mat_sz_1; i++) //must be square matrix
    {
      ierr=MatGetRow(mat,i,NULL,NULL,&mat_row);  CHKERRV(ierr);
      // use only the required part of the array  
      for(int j=0; j<vec_sz;j++)
      {
        mat_row_sbv[j]=mat_row[j];
        //LOG(DEBUG) << "mat_row[" << j << "]" << mat_row[j];
      }
           
      ierr=VecSetValues(mat_row_vec,vec_sz,idx,mat_row_sbv,INSERT_VALUES); CHKERRV(ierr);
      VecAssemblyBegin(mat_row_vec);
      VecAssemblyEnd(mat_row_vec);
      
      //VecGetSize(mat_row_vec,&vec_sz_tmp);
      //LOG(DEBUG) << "vec_sz_tmp: " << vec_sz_tmp;
      ierr=VecTDot(x,mat_row_vec,&val[i]);  CHKERRV(ierr);
      //LOG(DEBUG) << "val[" << i << "]: " << val[i];
    }
    
    ierr=VecSetValues(y,mat_sz_1,idx_2,val,INSERT_VALUES); CHKERRV(ierr); // would it work for parallel!?
    VecAssemblyBegin(y);
    VecAssemblyEnd(y);
    
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
    PetscInt *idx;
    PetscMalloc1(vec_sz,&idx);
    for(PetscInt i=0; i<vec_sz; i++)
      idx[i]=i;
    
    PetscInt *idx_2;
    PetscMalloc1(mat_sz_1,&idx_2);
    for(PetscInt i=0; i<mat_sz_2; i++) //square matrix
      idx_2[i]=i;
    
    const PetscScalar *mat_row;
    Vec mat_row_vec; //to allocate
    ierr=VecDuplicate(x,&mat_row_vec); CHKERRV(ierr);
    
    PetscScalar val[vec_sz];
    
    for(int i=0; i<vec_sz; i++)
    {
      ierr=MatGetRow(mat,i,NULL,NULL,&mat_row);  CHKERRV(ierr);      
      ierr=VecSetValues(mat_row_vec,mat_sz_2,idx_2,mat_row,INSERT_VALUES); CHKERRV(ierr);
      VecAssemblyBegin(mat_row_vec);
      VecAssemblyEnd(mat_row_vec);
      
      ierr=VecTDot(x,mat_row_vec,&val[i]);  CHKERRV(ierr);
    }
    
    ierr=VecSetValues(y,vec_sz,idx,val,INSERT_VALUES); CHKERRV(ierr); // would it work for parallel!?
    VecAssemblyBegin(y);
    VecAssemblyEnd(y);
  }
}

} //namespace
