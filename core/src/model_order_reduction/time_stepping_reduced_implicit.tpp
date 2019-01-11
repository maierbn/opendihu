#include "model_order_reduction/time_stepping_reduced_implicit.h"

#include <Python.h>
#include "utility/python_utility.h"
#include <petscksp.h>
#include "solver/solver_manager.h"
#include "solver/linear.h"
#include "data_management/time_stepping/time_stepping_implicit.h"

namespace ModelOrderReduction
{
template<typename TimeSteppingImplicitType>
TimeSteppingSchemeOdeReducedImplicit<TimeSteppingImplicitType>::
TimeSteppingSchemeOdeReducedImplicit(DihuContext context):
TimeSteppingSchemeOdeReduced<TimesteppingImplicitType>(context,"ImplicitEulerReduced"),
initialized_(false)
{
}

template<typename TimeSteppingImplicitType>
void TimeSteppingSchemeOdeReducedImplicit<TimeSteppingImplicitType>::
initialize()
{
  if (this->initialized_)
    return;
  
  LOG(TRACE) << "TimeSteppingSchemeOdeReducedImplicit::initialize()";
  
  // TO BE IMPLEMENTED
  TimeSteppingSchemeOdeReduced::initialize();
  
  // compute the system matrix
  this->setSystemMatrix(this->timeStepWidth_);
  
  // initialize the linear solver that is used for solving the implicit system
  initializeLinearSolver();
  
  // set matrix used for linear system and preconditioner to ksp context
  Mat &redSysMatrix = this->dataMOR_->redSysMatrix()->valuesGlobal();
  assert(this->ksp_);
  PetscErrorCode ierr;
  ierr = KSPSetOperators(*ksp_, redSysMatrix, redSysMatrix); CHKERRV(ierr);
  
  // TO BE IMPLEMENTED
  
  TimeSteppingSchemeOdeReduced<TimeSteppingImplicitType>::initialize();
  
  this->initialized_=true;
  
}

template<typename TimeSteppingImplicitType>
void ImplicitEulerReduced<TimeSteppingImplicitType>::
setSystemMatrix(double timeStepWidth)
{
  this->fullTimestepping_.setSystemMatrix(timeStepWidth);
  this->setRedSystemMatrix(this->fullTimestepping_.dataImplicit_->systemMatrix(), this->dataMOR_->redSystemMatrix());
}

template<typename FunctionSpaceRowsType>
void MORBase<FunctionSpaceRowsType>::
setRedSystemMatrix()
{
  //to be implemented
  Mat &basis = this->dataMOR_->basis()->valuesGlobal();
  Mat &basisTransp = this->dataMOR_->basisTransp()->valuesGlobal(); 
  
  std::shared_ptr<PartitionedPetscMat<typename this->fullTimestepping_.FunctionSpace>> ptr_systemMatrix=this->fullTimestepping_.data().systemMatrix();
  Mat systemMatrix=ptr_systemMatrix->valuesGlobal();
  Mat redSystemMatrix=this->dataMOR_.redSystemMatrix()->valuesGlobal();
  
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

template<typename TimeSteppingImplicitType>
void TimeSteppingSchemeOdeReducedImplicit<TimeSteppingImplicitType>::
solveLinearSystem(Vec &input, Vec &output)
{
  // solve systemMatrix*output = input for output
  Mat &redSysMatrix = this->dataMOR_->redSysMatrix()->valuesGlobal();
  
  PetscErrorCode ierr;
  PetscUtility::checkDimensionsMatrixVector(redSysMatrix, input);
  
  // solve the system, KSPSolve(ksp,b,x)
  ierr = KSPSolve(*ksp_, input, output); CHKERRV(ierr);
  
  int numberOfIterations = 0;
  PetscReal residualNorm = 0.0;
  ierr = KSPGetIterationNumber(*ksp_, &numberOfIterations); CHKERRV(ierr);
  ierr = KSPGetResidualNorm(*ksp_, &residualNorm); CHKERRV(ierr);
  
  KSPConvergedReason convergedReason;
  ierr = KSPGetConvergedReason(*ksp_, &convergedReason); CHKERRV(ierr);
  
  VLOG(1) << "Linear system of implicit time stepping solved in " << numberOfIterations << " iterations, residual norm " << residualNorm
  << ": " << PetscUtility::getStringLinearConvergedReason(convergedReason);
}

template<typename TimeSteppingImplicitType>
void TimeSteppingSchemeOdeReducedImplicit<TimeSteppingImplicitType>::
initializeLinearSolver()
{ 
  if (linearSolver_ == nullptr)
  {
    LOG(DEBUG) << "Implicit time stepping: initialize linearSolver";
    
    // retrieve linear solver
    linearSolver_ = this->context_.solverManager()->template solver<Solver::Linear>(
      this->specificSettings_, this->data_->functionSpace()->meshPartition()->mpiCommunicator());
    ksp_ = linearSolver_->ksp();
  }
  else 
  {
    VLOG(2) << ": linearSolver_ already set";
  }
}


  
} //namespace
