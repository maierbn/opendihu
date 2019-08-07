#include "model_order_reduction/model_order_reduction.h"
#include "data_management/data.h"
//#include <petscmat.h>
#include <array>
#include "utility/svd_utility.h"

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
    //ierr=MatShift(basisTransp, 1); CHKERRV(ierr); //identity matrix to check
    
    Mat &basis=this->dataMOR_->basis()->valuesGlobal();
    
    //ierr=MatShift(basis, 1); CHKERRV(ierr); //identity matrix to check
    
    PetscInt mat_sz_1, mat_sz_2;
    MatGetSize(basis,&mat_sz_1,&mat_sz_2);
    LOG(DEBUG) << "basis, mat_sz_1: " << mat_sz_1 << "basis, mat_sz_2: " << mat_sz_2 << "==============";
    
    // input data is the transpose of the snapshot matrix    
    std::string inputData = specificSettingsMOR_.getOptionString("snapshots","");
    cout << inputData;
    std::vector<double> parsedCSV = SvdUtility::readCSV(inputData);  
    
    // input data is the transpose of the snapshot matrix    
    int columnsSnapshots = SvdUtility::getCSVRowCount(inputData);
    if(columnsSnapshots != mat_sz_2)
      LOG(ERROR) << "There exists " << columnsSnapshots << " snapshots while there are " << mat_sz_2 << " modes required for the basis"; 
    
    // input data is the transpose of the snapshot matrix    
    int rowsSnapshots = SvdUtility::getCSVColumnCount(inputData);        
    if(rowsSnapshots < mat_sz_1)
      LOG(ERROR) << "snapshots have the length (spatial resolution) " << rowsSnapshots << " but the basis has the length " << mat_sz_1;
    
    // svd decomposition
    //std::vector<double> result = SvdUtility::getSVD(parsedCSV, rowsSnapshots, columnsSnapshots);
    /*
    for(std::vector<double>::iterator it = result.begin(); it!=result.end(); ++it)
    {    
      std::cout << ' ' << *it;
    }
    */
    
    double* v = new double[rowsSnapshots*columnsSnapshots];
    double* v_reconst = new double[rowsSnapshots*columnsSnapshots];
    
    std::copy(parsedCSV.begin(), parsedCSV.end(), v); // This is the transpose of the snapshot matrix    
    //int cnt=0;
    //for(std::vector<double>::iterator it = parsedCSV.begin(); it!=parsedCSV.end(); ++it)
    //{    
    //  std::cout << ' ' << *it;
    //  v[cnt]=*it;
    //  cnt++;
    //}
    //std::cout << std::endl;
    
    int n = std::min(rowsSnapshots, columnsSnapshots);
    
    // J > K => only the first K cols of U are computed
    double* leftSingVec = new double[rowsSnapshots * n];
    
    // J < K => only the first J rows of T^T are computed
    double* rightSingVecT = new double[n * columnsSnapshots];
    
    //double* sigmas = new double[n];
    double* sigma = new double[n * n];
    
    // snapshots are required in the column-major-format, therefore, v is the transpose of the snapshots as csv data
    SvdUtility::getSVD(v,rowsSnapshots,columnsSnapshots,leftSingVec,sigma,rightSingVecT);
    SvdUtility::reconstructSnapshots(rowsSnapshots, columnsSnapshots, leftSingVec, sigma,rightSingVecT, v_reconst);
    
    PetscInt *idx_1;
    PetscMalloc1(mat_sz_1,&idx_1);
    for(PetscInt i=0; i<mat_sz_1; i++)
      idx_1[i]=i;
    
    PetscInt *idx_2;
    PetscMalloc1(mat_sz_2,&idx_2);
    for(PetscInt i=0; i<mat_sz_2; i++)
      idx_2[i]=i;
    
    ierr=MatSetValues(basis,mat_sz_1,idx_1,mat_sz_2,idx_2,leftSingVec,INSERT_VALUES); CHKERRV(ierr);
    MatAssemblyBegin(basis,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(basis,MAT_FINAL_ASSEMBLY);
    
    PetscScalar val;
    
    cout << "basis" << endl;
    for (PetscInt i=0; i< mat_sz_1; i++)
    {
      for(PetscInt j=0; j<mat_sz_2; j++)
      {
        MatGetValues(basis,1,&i,1,&j,&val);
        MatSetValues(basisTransp,1,&j,1,&i,&val,INSERT_VALUES);
        cout << val;
      }
      cout << std::endl;
    }
    MatAssemblyBegin(basisTransp,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(basisTransp,MAT_FINAL_ASSEMBLY);
    
    cout << endl;
    
    //ierr=MatTranspose(basis,MAT_INITIAL_MATRIX,&basisTransp); CHKERRV(ierr); //produces err
    
    //PetscScalar val;
    cout << "basisTransp" << endl;
    for (PetscInt i=0; i< mat_sz_2; i++)
    {
      for(PetscInt j=0; j<mat_sz_1; j++)
      {
        MatGetValues(basisTransp,1,&i,1,&j,&val);
        cout << val;
      }
      cout << std::endl;
    }
    /*
    PetscInt mat_sz_1, mat_sz_2;
    MatGetSize(basisTransp,&mat_sz_1,&mat_sz_2);
    LOG(DEBUG) << "mat_sz_1: " << mat_sz_1 << "mat_sz_2" << mat_sz_2 << "==============";
    */       
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
    PetscInt vec_sz,vec_sz_y,mat_sz_1,mat_sz_2;
    
    VecGetSize(x,&vec_sz); // size of the full-order solution
    VecGetSize(y,&vec_sz_y); // size of the reduced solution
    MatGetSize(mat,&mat_sz_1,&mat_sz_2);
    
    VLOG(2) << "mat_sz_2: " << mat_sz_2 << " vec_sz: " << vec_sz;
    VLOG(2) << "mat_sz_1: " << mat_sz_1 << " vec_sz_y: " << vec_sz_y;
    
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
          ierr=VecSetValues(mat_row_vec,vec_sz,idx,mat_row,INSERT_VALUES); CHKERRV(ierr);
          VecAssemblyBegin(mat_row_vec);
          VecAssemblyEnd(mat_row_vec);
          
          ierr=VecTDot(mat_row_vec,x,&val[i]);  CHKERRV(ierr);
        }
        
        ierr=VecSetValues(y,mat_sz_1,idx_2,val,INSERT_VALUES); CHKERRV(ierr); // would it work for parallel!?
        VecAssemblyBegin(y);
        VecAssemblyEnd(y);    
    }
    else
      LOG(ERROR) << "smaller size of the out put vector " << vec_sz << "than matrix size " <<  mat_sz_2 << " in matrix vector multiplication.";
    }
  }

  template<typename FunctionSpaceRowsType>
  void MORBase<FunctionSpaceRowsType>::
  MatMultFull(Mat mat,Vec x,Vec y)
  {
    PetscErrorCode ierr;  
    PetscInt vec_sz_x,vec_sz_y,mat_sz_1,mat_sz_2;
    
    VecGetSize(x,&vec_sz_x); // size of the full-order solution
    VecGetSize(y,&vec_sz_y); // size of the full-order solution
    MatGetSize(mat,&mat_sz_1,&mat_sz_2);
    
    VLOG(2) << "mat_sz_2: " << mat_sz_2 << " vec_sz_x: " << vec_sz_x;
    VLOG(2) << "mat_sz_1: " << mat_sz_2 << " vec_sz_y: " << vec_sz_y;
    
    // if the full-order and reduced-order solutions have the same length!?
    // This is not the case for the diffusion term in electrophysiology. 
    if(mat_sz_1==vec_sz_y) 
    {
      if(mat_sz_2==vec_sz_x) // compatibility of the multiplication
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
      PetscMalloc1(vec_sz_y,&idx);
      for(PetscInt i=0; i<vec_sz_y; i++)
        idx[i]=i;
      
      PetscInt *idx_2;
      PetscMalloc1(mat_sz_2,&idx_2);
      for(PetscInt i=0; i<mat_sz_2; i++)
        idx_2[i]=i;
      
      const PetscScalar *mat_row;
      Vec mat_row_vec;
      ierr=VecDuplicate(x,&mat_row_vec); CHKERRV(ierr);
      
      PetscScalar *val;
      PetscMalloc1(vec_sz_y,&val);
      
      for(int i=0; i<vec_sz_y; i++) // not all rows of the mat (basis) would be multiplied because size of y (full-order solution) could be smaller than rows of mat.
      {
        ierr=MatGetRow(mat,i,NULL,NULL,&mat_row);  CHKERRV(ierr);      
        ierr=VecSetValues(mat_row_vec,mat_sz_2,idx_2,mat_row,INSERT_VALUES); CHKERRV(ierr);
        VecAssemblyBegin(mat_row_vec);
        VecAssemblyEnd(mat_row_vec);
        
        ierr=VecTDot(mat_row_vec,x,&val[i]);  CHKERRV(ierr);
      }
      
      ierr=VecSetValues(y,vec_sz_y,idx,val,INSERT_VALUES); CHKERRV(ierr); // would it work for parallel!?
      VecAssemblyBegin(y);
      VecAssemblyEnd(y);
    }
}

} //namespace
