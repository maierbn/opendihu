#include "spatial_discretization/finite_element_method/05_time_stepping.h"

#include <Python.h>
#include <iostream>
#include <petscmat.h>
#include <petscksp.h>
#include <vector>
#include <numeric>

#include "easylogging++.h"

#include "control/types.h"
#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "solver/solver_manager.h"
#include "solver/linear.h"


namespace SpatialDiscretization
{

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
setInvLumMassMatrix()
{
  const int D = BasisOnMeshType::dim();
  LOG(TRACE)<<"setInvLumMassMatrix" << D << "D";
  
  if (!this->data_.invLumMassMatrixInitialized())
      this->data_.initializeInvLumMassMatrix();
  
  Mat &invLumMassMatrix=this->data_.invLumMassMatrix();
  Mat &massMatrix=this->data_.massMatrix();
  
  PetscErrorCode ierr;
  
  PetscInt n_row, n_col;
  ierr=MatGetSize(massMatrix,&n_row,&n_col);
  LOG(INFO)<<"n_row "<< n_row <<" n_col " << n_col;
     
  Vec rowSum; 
  ierr=VecCreate(PETSC_COMM_WORLD,&rowSum);
  ierr=VecSetType(rowSum,VECMPI);  
  ierr=VecSetSizes(rowSum,PETSC_DECIDE,n_row);
  ierr=VecAssemblyBegin(rowSum);
  ierr=VecAssemblyEnd(rowSum);
 
  //In case of linear and bilinear basis functions
  ierr=MatGetRowSum(massMatrix,rowSum);
  //for the inverse matrix
  ierr=VecReciprocal(rowSum); 
 
  /*
  PetscInt rowSum_size;
  ierr=VecGetSize(rowSum,&rowSum_size);
  LOG(INFO)<<"rowSum_size"<< rowSum_size;
  
  PetscScalar *rowSum_value[rowSum_size];
  ierr=VecGetArray(rowSum,rowSum_value);
  */

  LOG(TRACE)<<"StartAssembleInvLumMassMatrix" << D << "D";

  std::shared_ptr<BasisOnMeshType> mesh = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh());
  dof_no_t n = mesh->nDofs();
  LOG(INFO)<< "dof_no_t "<< n;
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  
  // set values
  int cntr = 0;
  // loop over elements
  for (element_no_t elementNo = 0; elementNo < mesh->nElements(); elementNo++)
  {
    auto dof = mesh->getElementDofNos(elementNo);

    for (int i=0; i<nDofsPerElement; i++)
    {
      for (int j=0; j<nDofsPerElement; j++)
      {
          VLOG(3) << " inversed lumped massMatrix entry ( " << dof[i] << "," << dof[j] << ") (no. " << cntr++ << ")";
          //LOG(DEBUG) << " inversed lumped massMatrix entry ( " << dof[i] << "," << dof[j] << ") (no. " << cntr++ << ")";
          ierr = MatSetValue(invLumMassMatrix, dof[i], dof[j], 0, INSERT_VALUES); CHKERRV(ierr);
      }
            
    }     
  }    
  
  double val_get;
  for (int i=0; i<n_row; i++)
  {
    //ierr=MatSetValue(invLumMassMatrix, i,i, *rowSum_value[i], INSERT_VALUES);
    
    ierr=VecGetValues(rowSum,1,&i,&val_get);
    //LOG(INFO)<<"val_get "<< val_get;   
    ierr=MatSetValue(invLumMassMatrix, i,i, val_get, INSERT_VALUES);
  }
  
    
  //For finite element basis greater than 2, to be implemented here
    
  ierr = MatAssemblyBegin(invLumMassMatrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(invLumMassMatrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);        
    
  this->invLumMassMatrixSet_=true;
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
bool FiniteElementMethodTimeStepping<BasisOnMeshType,QuadratureType,Term>::
invLumMassMatrixSet(){
  return this->invLumMassMatrixSet_;
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
preComputeSystemMatrix(double timeStepWidth)
{
  const int D = BasisOnMeshType::dim();
  LOG(TRACE)<<"preComputeSystemMatrix" << D << "D";
  
  if(!this->invLumMassMatrixSet())
    this->setInvLumMassMatrix();
 
  Mat &systemMatrix = this->data_.systemMatrix();
  Mat &invLumMassMatrix = this->data_.invLumMassMatrix();
  Mat &stiffnessMatrix = this->data_.stiffnessMatrix();
  
  const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  
  PetscErrorCode ierr;
  std::shared_ptr<BasisOnMeshType> mesh = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh());
  
  ierr=MatMatMult(invLumMassMatrix, stiffnessMatrix,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&systemMatrix);
  PetscScalar scale=-timeStepWidth;
  ierr=MatScale(systemMatrix,scale);
  
  // set diagonal values
  int cntr = 1;
  // loop over elements
  for (element_no_t elementNo = 0; elementNo < mesh->nElements(); elementNo++)
  {
    auto dof = mesh->getElementDofNos(elementNo);
    
    for (int i=0; i<nDofsPerElement; i++)
    {
      VLOG(3) << " systemMatrix entry ( " << dof[i] << "," << dof[i] << ") (no. " << cntr++ << ")";
      //LOG(DEBUG) << " systemMatrix entry ( " << dof[i] << "," << dof[i] << ") (no. " << cntr++ << ")";
      ierr = MatSetValue(systemMatrix, dof[i], dof[i], 1.0, ADD_VALUES); CHKERRV(ierr);      
    }
  }
  
  ierr = MatAssemblyBegin(systemMatrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(systemMatrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);  
  
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
solveLinearSystem(Vec &input, Vec &output)
{
  //Au^(t+1)=u^(t)
  
  Mat &systemMatrix=this->data_.systemMatrix();
  
  PetscErrorCode ierr;
  
  // create linear solver context
  std::shared_ptr<Solver::Linear> linearSolver = this->context_.solverManager()->template solver<Solver::Linear>(this->specificSettings_);
  std::shared_ptr<KSP> ksp = linearSolver->ksp();
  
  // set matrix used for linear system and preconditioner to ksp context
  ierr = KSPSetOperators (*ksp, systemMatrix, systemMatrix); CHKERRV(ierr);
  
  //use the default initial guess (zero) by the PETSC because the input and output are the same by the implicit method   
  // solve the system
  ierr = KSPSolve(*ksp, input, output); CHKERRV(ierr);
  
  int numberOfIterations = 0;
  PetscReal residualNorm = 0.0;
  ierr = KSPGetIterationNumber(*ksp, &numberOfIterations); CHKERRV(ierr);
  ierr = KSPGetResidualNorm(*ksp, &residualNorm); CHKERRV(ierr);

  //LOG(INFO) << "Rhs recovered in " << numberOfIterations << " iterations, residual norm " << residualNorm;
  VLOG(1) << "Solution of the linear system recovered in " << numberOfIterations << " iterations, residual norm " << residualNorm;
  
}

} // namespace SpatialDiscretization