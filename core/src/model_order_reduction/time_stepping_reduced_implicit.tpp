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
  
  // TO BE IMPLEMENTED
  
  TimeSteppingSchemeOdeReduced<TimeSteppingImplicitType>::initialize();
  
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
  //to be implemented
  Mat &basis = this->dataMOR_->basis()->valuesGlobal();
  Mat &basisTransp = this->dataMOR_->basisTransp()->valuesGlobal(); 
   
  Mat &systemMatrix=this->fullTimestepping_.dataImplicit().systemMatrix()->valuesGlobal();
  Mat &redSystemMatrix=this->dataMOR_->redSystemMatrix()->valuesGlobal();
  
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
    std::shared_ptr<Partition::MeshPartition<typename TimeSteppingImplicitType::FunctionSpace>> 
    meshPartitionRows=this->fullTimestepping_.dataImplicit().systemMatrix()->meshPartitionRows(); //compatible for multiplication  
    
    std::shared_ptr<Partition::MeshPartition<typename ::FunctionSpace::Generic>>
    meshPartitionColumns=this->dataMOR_->basisTransp()->meshPartitionRows();    
    
    //std::shared_ptr<PartitionedPetscMat<typename ::FunctionSpace::Generic,typename TimeSteppingImplicitType::FunctionSpace>> basis_submatrix;
    std::shared_ptr<PartitionedPetscMat<typename ::FunctionSpace::Generic,typename TimeSteppingImplicitType::FunctionSpace>> basisTransp_submatrix;
    
    //basis_submatrix = std::make_shared<PartitionedPetscMat<::FunctionSpace::Generic,typename TimeSteppingImplicitType::FunctionSpace>>(
      //meshPartitionRows, meshPartitionColumns, 1, "basis_submatrix");
    basisTransp_submatrix= std::make_shared<PartitionedPetscMat<::FunctionSpace::Generic,typename TimeSteppingImplicitType::FunctionSpace>>(
      meshPartitionColumns, meshPartitionRows, 1, "basisTransp_submatrix");
      
    //Mat &basis_sbm=basis_submatrix->valuesGlobal();
    Mat basis_sbm;
    Mat basisTransp_sbm=basisTransp_submatrix->valuesGlobal();
    //Mat &basisTransp_sbm;
    
    const PetscScalar *vals_total;
    //const PetscScalar vals_total[ncols_bs];
    //const PetscScalar vals_sbm[nrows_sys];
    
    for( int row=0; row<nrows_bs;row++)
    {
      // to get each row of the V_k^T
      MatGetRow(basisTransp,row,NULL,NULL,&vals_total); CHKERRV(ierr);    
      MatSetValuesRow(basisTransp_sbm,row,vals_total); CHKERRV(ierr); //inconsistent sizes may work. Otherwise build vals_sbm.     
    }
    
    MatCreateTranspose(basisTransp_sbm,&basis_sbm); CHKERRV(ierr);
    
    //D=A*B*C
    MatMatMatMult(basisTransp_sbm,systemMatrix,basis_sbm,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&redSystemMatrix); CHKERRV(ierr);
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
