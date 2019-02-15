#include "model_order_reduction/mor.h"
#include "data_management/data.h"
//#include <petscmat.h>
#include <array>

namespace ModelOrderReduction
{

  template<typename FunctionSpaceRowsType>
  MORBase<FunctionSpaceRowsType>::
  MORBase(DihuContext context):
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

    // assemble matrices for parallel use, see https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatAssemblyBegin.html
    this->dataMOR_->basisTransp()->assembly(MAT_FINAL_ASSEMBLY);
    this->dataMOR_->basis()->assembly(MAT_FINAL_ASSEMBLY);

    Mat &basisTransp=this->dataMOR_->basisTransp()->valuesGlobal();
    
    PetscErrorCode ierr;
    ierr=MatShift(basisTransp, 1); CHKERRV(ierr); //identity matrix to check
    
    Mat &basis=this->dataMOR_->basis()->valuesGlobal();
    
    ierr=MatShift(basis, 1); CHKERRV(ierr); //identity matrix to check
    
    
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

  template<typename FunctionSpaceRowsType>
  void MORBase<FunctionSpaceRowsType>::
  MatMultReduced(Mat mat,Vec x,Vec y)
  {
    PetscErrorCode ierr;
    PetscInt vec_sz,mat_sz_1,mat_sz_2;
    
    VecGetSize(x,&vec_sz); // size of the full-order solution
    MatGetSize(mat,&mat_sz_1,&mat_sz_2);
    
    //LOG(TRACE) << "mat_sz_2: " << mat_sz_2 << " vec_sz: " << vec_sz;
    
    if(mat_sz_2==vec_sz)
    {
      ierr=MatMult(mat,x,y); CHKERRV(ierr);
    }
    else
    {
      if(mat_sz_2 > vec_sz) // mat_sz_2 includes all components must be greater than only one variable (transmembrane potential)
      {
        //LOG(TRACE) << "MatMultReduced";
        
        // columns of mat (transpose of the basis vector) are in fact equal to the 
        // rows of the snapshot matrix (containing all components in total 
        // reduction of electrophysiology), which could be more than the size of 
        // the full-order solution vector of the transmembrane potential for the
        // diffusion term in electrophysiology. mat_sz_2 is equal to 4n 
        // (Hodgkin-Huxley) but vec_sz is n (size of the action potential). 
        
        PetscInt *idx;
        PetscMalloc1(vec_sz,&idx);
        for(PetscInt i=0; i<vec_sz; i++)
              idx[i]=i;
        
        PetscInt *idx_2;
        PetscMalloc1(mat_sz_1,&idx_2);
        for(PetscInt i=0; i<mat_sz_1; i++)
          idx_2[i]=i;
              
        const PetscScalar *mat_row;
        
        Vec mat_row_vec;
        ierr=VecDuplicate(x,&mat_row_vec); CHKERRV(ierr);
        
        PetscScalar *val;
        PetscMalloc1(mat_sz_1,&val);
          
        for(int i=0; i<mat_sz_1; i++)
        {
          ierr=MatGetRow(mat,i,NULL,NULL,&mat_row);  CHKERRV(ierr); // can take only the rows of the current processor
          ierr=VecSetValuesLocal(mat_row_vec,vec_sz,idx,mat_row,INSERT_VALUES); CHKERRV(ierr);
          VecAssemblyBegin(mat_row_vec);
          VecAssemblyEnd(mat_row_vec);
          
          ierr=VecTDot(mat_row_vec,x,&val[i]);  CHKERRV(ierr);
        }
        
        ierr=VecSetValues(y,mat_sz_1,idx_2,val,INSERT_VALUES); CHKERRV(ierr); // would it work for parallel!?
        VecAssemblyBegin(y);
        VecAssemblyEnd(y);    
    }
    else
      LOG(ERROR) << "smaller size of the out put vector in matrix vector multiplication.";
    }
  }

  template<typename FunctionSpaceRowsType>
  void MORBase<FunctionSpaceRowsType>::
  MatMultFull(Mat mat,Vec x,Vec y)
  {
    PetscErrorCode ierr;  
    PetscInt vec_sz,mat_sz_1,mat_sz_2;
    
    VecGetSize(y,&vec_sz); // size of the full-order solution
    MatGetSize(mat,&mat_sz_1,&mat_sz_2);
    
    //LOG(TRACE) << "mat_sz_2: " << mat_sz_2 << " vec_sz: " << vec_sz;
    
    // if the full-order and reduced-order solutions have the same length!?
    // This is not the case for the diffusion term in electrophysiology. 
    if(mat_sz_1==vec_sz) 
    {
      if(mat_sz_2==vec_sz) // compatibility of the multiplication
      {
        ierr=MatMult(mat,x,y); CHKERRV(ierr);
      }
      else 
        LOG(ERROR) << "incompatible size for matrix vector multiplication";
    }
    else
    {  
      //LOG(TRACE) << "MatMultFull";
      
      PetscInt *idx;
      PetscMalloc1(vec_sz,&idx);
      for(PetscInt i=0; i<vec_sz; i++)
        idx[i]=i;
      
      PetscInt *idx_2;
      PetscMalloc1(mat_sz_2,&idx_2);
      for(PetscInt i=0; i<mat_sz_2; i++)
        idx_2[i]=i;
      
      const PetscScalar *mat_row;
      Vec mat_row_vec;
      ierr=VecDuplicate(x,&mat_row_vec); CHKERRV(ierr);
      
      PetscScalar *val;
      PetscMalloc1(vec_sz,&val);
      
      for(int i=0; i<vec_sz; i++) // not all rows of the mat (basis) would be multiplied because size of y (full-order solution) could be smaller than rows of mat.
      {
        ierr=MatGetRow(mat,i,NULL,NULL,&mat_row);  CHKERRV(ierr);      
        ierr=VecSetValuesLocal(mat_row_vec,mat_sz_2,idx_2,mat_row,INSERT_VALUES); CHKERRV(ierr);
        VecAssemblyBegin(mat_row_vec);
        VecAssemblyEnd(mat_row_vec);
        
        ierr=VecTDot(mat_row_vec,x,&val[i]);  CHKERRV(ierr);
      }
      
      ierr=VecSetValuesLocal(y,vec_sz,idx,val,INSERT_VALUES); CHKERRV(ierr); // would it work for parallel!?
      VecAssemblyBegin(y);
      VecAssemblyEnd(y);
    }
}

} //namespace
