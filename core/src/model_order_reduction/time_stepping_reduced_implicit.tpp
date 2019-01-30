#include "model_order_reduction/time_stepping_reduced_implicit.h"

#include <Python.h>
#include "utility/python_utility.h"
#include <petscksp.h>
#include <petscmat.h> 
#include "solver/solver_manager.h"
#include "solver/linear.h"
#include "data_management/time_stepping/time_stepping_implicit.h"
#include "time_stepping_scheme/time_stepping_scheme_ode.h"

namespace ModelOrderReduction
{
template<typename TimeSteppingImplicitType>
TimeSteppingSchemeOdeReducedImplicit<TimeSteppingImplicitType>::
TimeSteppingSchemeOdeReducedImplicit(DihuContext context,std::string name):
TimeSteppingSchemeOdeReduced<TimeSteppingImplicitType>(context,name)
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
  TimeSteppingSchemeOdeReduced<TimeSteppingImplicitType>::initialize();
  
  // compute the system matrix
  this->setSystemMatrix(this->timeStepWidth_);
  
  // set the boundary conditions to system matrix, i.e. zero rows and columns of Dirichlet BC dofs and set diagonal to 1
  this->fullTimestepping_.dirichletBoundaryConditions()->applyInSystemMatrix(this->fullTimestepping_.dataImplicit().systemMatrix(), this->fullTimestepping_.dataImplicit().boundaryConditionsRightHandSideSummand());
  
  // initialize the linear solver that is used for solving the implicit system
  initializeLinearSolver();
  
  // set matrix used for linear system and preconditioner to ksp context
  Mat redSystemMatrix = this->dataMOR_->redSystemMatrix()->valuesGlobal();
  assert(this->ksp_);
  PetscErrorCode ierr;
  ierr = KSPSetOperators(*ksp_, redSystemMatrix, redSystemMatrix); CHKERRV(ierr);
  
  this->initialized_=true;
  
}

template<typename TimeSteppingImplicitType>
void TimeSteppingSchemeOdeReducedImplicit<TimeSteppingImplicitType>::
setSystemMatrix(double timeStepWidth)
{
  // is it set by initializing the fullTimestepping?
  //this->fullTimestepping_.setSystemMatrix(timeStepWidth);
  this->setRedSystemMatrix();
}

template<typename TimeSteppingImplicitType>
void TimeSteppingSchemeOdeReducedImplicit<TimeSteppingImplicitType>::
setRedSystemMatrix()
{

  Mat &basis = this->dataMOR_->basis()->valuesGlobal();
  Mat &basisTransp = this->dataMOR_->basisTransp()->valuesGlobal(); 
  
  Mat &systemMatrix=this->fullTimestepping_.dataImplicit().systemMatrix()->valuesGlobal();
  Mat &redSystemMatrix=this->dataMOR_->redSystemMatrix()->valuesGlobal();
  
  PetscErrorCode ierr; 
    
  PetscInt mat_sz_1, mat_sz_2;
  
  MatGetSize(basisTransp,&mat_sz_1,&mat_sz_2);
  PetscInt nrows_bs=mat_sz_1;
  PetscInt ncols_bs=mat_sz_2;
  LOG(DEBUG) << "setRedSystemMatrix, nrows_bs: " << nrows_bs << " ncols_bs: " << ncols_bs;
  
  MatGetSize(systemMatrix,&mat_sz_1,&mat_sz_2);
  PetscInt nrows_sys=mat_sz_1; //square matrix
  PetscInt ncols_sys=mat_sz_2;
  
  LOG(DEBUG) << "setRedSystemMatrix, nrows_sys: " << nrows_sys << " ncols_sys: " << ncols_sys;
  
  if(ncols_bs == nrows_sys)
  {
    //! Reduction of the system matrix in case of compatible row spaces of the system matrix and reduced basis    
    //D=A*B*C
    //This method is not able to multiply combination of three sparse and dense matrices
    //MatMatMatMult(basisTransp,systemMatrix,basis,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&redSystemMatrix); CHKERRV(ierr);
    Mat matrix;
    MatDuplicate(basis,MAT_DO_NOT_COPY_VALUES,&matrix);
    MatMatMult(systemMatrix,basis,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&matrix); CHKERRV(ierr);
    MatMatMult(basisTransp,matrix,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&redSystemMatrix); CHKERRV(ierr);
  }
  else
  {    
    std::shared_ptr<Partition::MeshPartition<typename TimeSteppingImplicitType::FunctionSpace>> 
    meshPartitionRows=this->dataMOR_->basis()->meshPartitionRows();
    
    std::shared_ptr<Partition::MeshPartition<typename ::FunctionSpace::Generic>>
    meshPartitionColumns=this->dataMOR_->basisTransp()->meshPartitionColumns();    
    
    std::shared_ptr<PartitionedPetscMat<::FunctionSpace::Generic,typename TimeSteppingImplicitType::FunctionSpace>> systemMatrix_ext;
    
    systemMatrix_ext= std::make_shared<PartitionedPetscMat<::FunctionSpace::Generic,::FunctionSpace::Generic>>(
      meshPartitionRows, meshPartitionColumns, 1, "systemMatrix_ext");          
    
    Mat &matrix_ext=systemMatrix_ext->valuesGlobal();    
    ierr=MatShift(matrix_ext,1); CHKERRV(ierr);
    
    PetscScalar *vals;
    PetscMalloc1(nrows_sys*nrows_sys,&vals);
    
    PetscInt *idx;
    PetscMalloc1(nrows_sys,&idx);
    for( int i=0; i<nrows_sys;i++)
    {
      idx[i]=i;
    }
    
    ierr=MatGetValues(systemMatrix,nrows_sys,idx,nrows_sys,idx,vals); CHKERRV(ierr);
    ierr=MatSetValuesLocal(matrix_ext,nrows_sys,idx,nrows_sys,idx,vals,INSERT_VALUES); CHKERRV(ierr); //would replace the extra ones on the diagonal
    
    MatAssemblyBegin(matrix_ext,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(matrix_ext,MAT_FINAL_ASSEMBLY);
    
    /*
    const PetscScalar *vals_ext; 
    for(int i=0; i<nrows_bs;i++)
    {
      MatGetRow(matrix_ext,i,NULL,NULL,&vals_ext);
      for(int j=0; j<nrows_bs; j++)
        LOG(DEBUG) << "matrix_ext[" << i << "," << j<< "]: " << vals_ext[j];
    } 
    */
    
    //D=A*B*C
    //This method is not able to multiply combination of three sparse and dense matrices
    //MatMatMatMult(basisTransp_sbm,systemMatrix,basis_sbm,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&redSystemMatrix); CHKERRV(ierr);
    Mat matrix;
    MatDuplicate(basis,MAT_DO_NOT_COPY_VALUES,&matrix);
    MatMatMult(matrix_ext,basis,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&matrix); CHKERRV(ierr);
        
    MatGetSize(basisTransp,&mat_sz_1,&mat_sz_2);
    LOG(DEBUG) << "basisTransp mat_sz_1, mat_sz_2: " << mat_sz_1 << " " << mat_sz_2;
    MatGetSize(matrix,&mat_sz_1,&mat_sz_2);
    LOG(DEBUG) << "matrix mat_sz_1, mat_sz_2: " << mat_sz_1 << " " << mat_sz_2;
    
    
    MatMatMult(basisTransp,matrix,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&redSystemMatrix); CHKERRV(ierr);
    
    /*
    const PetscScalar *vals_total;    
    for(int i=0; i<nrows_bs;i++)
    {
      MatGetRow(redSystemMatrix,i,NULL,NULL,&vals_total);
      for(int j=0; j<nrows_bs; j++)
        LOG(DEBUG) << "redSystemMatrix[" << i << "," << j<< "]: " << vals_total[j];
    }
    */
       
  }
  
}


template<typename TimeSteppingImplicitType>
void TimeSteppingSchemeOdeReducedImplicit<TimeSteppingImplicitType>::
solveLinearSystem(Vec &input, Vec &output)
{
  // solve systemMatrix*output = input for output
  Mat &redSystemMatrix = this->dataMOR_->redSystemMatrix()->valuesGlobal();
  
  PetscErrorCode ierr;
  PetscUtility::checkDimensionsMatrixVector(redSystemMatrix, input);
  
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
