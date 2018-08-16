#include "finite_element_method_time_stepping/05_finite_element_method_time_stepping.h"

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
evaluateTimesteppingRightHandSideImplicit(Vec &input, Vec &output, int timeStepNo, double currentTime)
{
  LOG(TRACE)<<"evaluateTimesteppingRightHandSideImplicit";
  
  // this method computes output = M*input  
  Mat &massMatrix=this->data_.massMatrix();
  
  PetscErrorCode ierr;
  
  /*
  PetscInt nEntries;
  VecGetSize(input, &nEntries);  
  PetscScalar v=0.0;  
  for (PetscInt i=0;i<nEntries;i++)
    for(PetscInt j=0;j<nEntries;j++)
    {      
      ierr=MatGetValues(massMatrix,1,&i,1,&j,&v);
      LOG(DEBUG)<<"val_get: massMatrix "<< v; 
    }
  */
  
  /*
  for (int i=0;i<nEntries;i++)
  {
    ierr=VecGetValues(input,1,&i,&val_get);
    LOG(DEBUG)<<"val_get: input "<< val_get;   
  }
  */
  
  checkDimensions(massMatrix,input);
  ierr=MatMult(massMatrix,input,output);
  
  /*
  for (int i=0;i<nEntries;i++)
  {
    ierr=VecGetValues(output,1,&i,&val_get);
    LOG(DEBUG)<<"val_get: output "<< val_get;   
  }
  */
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
preComputeSystemMatrix1()
{
  const int D = BasisOnMeshType::dim();
  LOG(TRACE)<<"preComputeSystemMatrix" << D << "D";  
  
  //modifying the stiffness matrix to use it as the system matrix
  this->data_.systemMatrix() = this->data_.stiffnessMatrix();
  Mat &massMatrix = this->data_.massMatrix();
  
  PetscErrorCode ierr;

  //Scale is 1.0 because the mass matrix is scaled with the time step by the initialization
  PetscScalar scale=1.0;
  ierr=MatAXPY(this->data_.systemMatrix(),scale,massMatrix,SAME_NONZERO_PATTERN);
  
}

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
    ierr=VecGetValues(rowSum,1,&i,&val_get);
    //LOG(DEBUG)<<"val_get "<< val_get;   
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
preComputeSystemMatrix(Mat &systemMatrix)
{
  const int D = BasisOnMeshType::dim();
  LOG(TRACE)<<"preComputeSystemMatrix" << D << "D";
  
  if(!this->invLumMassMatrixSet())
    this->setInvLumMassMatrix();
  
  Mat &invLumMassMatrix = this->data_.invLumMassMatrix();
  Mat &stiffnessMatrix = this->data_.stiffnessMatrix();
  
  PetscErrorCode ierr;
       
  //The result matrix is created by the routine itself
  ierr=MatMatMult(invLumMassMatrix, stiffnessMatrix,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&systemMatrix);
 
  
  PetscInt n_row, n_col;
  ierr=MatGetSize(stiffnessMatrix,&n_row,&n_col);
  LOG(DEBUG) << "n_row: " << n_row;
   
 /*
  //Scale is 1.0 because the mass matrix is scaled with the time step by the initialization
  PetscScalar scale=1.0;
  ierr=MatScale(systemMatrix,scale);
*/  
  
  for (int i=0; i<n_row; i++)
    ierr = MatSetValue(systemMatrix, i, i, 1.0, ADD_VALUES); CHKERRV(ierr);
  
  /*
   P ets*cScalar v=0.0;  
   for (PetscInt i=0;i<n_row;i++)
     for(PetscInt j=0;j<n_col;j++)
     {          
        ierr=MatGetValues(systemMatrix,1,&i,1,&j,&v);
        LOG(INFO)<<"val_get: systemMatrix "<< v;  
     }    
  */
  
  ierr = MatAssemblyBegin(systemMatrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(systemMatrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);    
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
solveLinearSystem(Vec &input, Vec &output)
{
  //Au^(t+1)=u^(t)
  
  Mat &systemMatrix=this->data_.systemMatrix();
  
  LOG(DEBUG) << "solveLinearSystem";
  
  PetscErrorCode ierr;
  
  /*
  PetscInt n_row, n_col;
  ierr=MatGetSize(systemMatrix,&n_row,&n_col);
  PetscScalar v=0.0; 
  for (PetscInt i=0;i<n_row;i++)
    for(PetscInt j=0;j<n_col;j++)
    {
      ierr=MatGetValues(systemMatrix,1,&i,1,&j,&v);
      LOG(DEBUG)<<"val_get: systemMatrix in solveLinearSystem: "<< v; 
    }
    */
   
   
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

  LOG(INFO) << "Linear system solved in " << numberOfIterations << " iterations, residual norm " << residualNorm;
  VLOG(1) << "Solution of the linear system recovered in " << numberOfIterations << " iterations, residual norm " << residualNorm;
  
}

} // namespace SpatialDiscretization