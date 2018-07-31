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
Mat &FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
setInvLumMassMatrix()
{
  
  // check if matrix discretization matrix exists
  if (!this->data_.invLumMassMatrixInitialized())
  {
    this->data_.initializeInvLumMassMatrix();
   
    const int D = BasisOnMeshType::dim();
    LOG(TRACE)<<"inverseLumpedMassMatrix" << D << "D";
  
  
    PetscErrorCode ierr;
  
    Mat &invLumMassMatrix=this->data_.inversedLumpedMassMatrix_;
    Mat &massMatrix=this->data_.massMatrix; 
    //for linear and bilinear basis functions
    Vec &rowSum;  
    ierr=MatGetRowSum(massMatrix,rowSum);
  
    PetscInt rowSum_size;
    ierr=VecGetSize(rowSum,rowSum_size);
  
    const int nDofsPerElement = BasisOnMeshType::nDofsPerElement();
  
    std::shared_ptr<BasisOnMeshType> mesh = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh());
  
    //for the inverse matrix
    ierr=VecReciprocal(rowSum);
  
    PetscScalar rowSum_value[rowSum_size];
    ierr=VecGetArray(rowSum,rowSum_value);

    // set values
    int cntr = 1;
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
        VLOG(3) << " lumped massMatrix entry ( " << dof[i] << "," << dof[i] << ") (no. " << cntr++ << ")";
        //LOG(DEBUG) << " inversed lumped massMatrix entry ( " << dof[i] << "," << dof[i] << ") (no. " << cntr++ << ")";
        ierr = MatSetValue(invLumMassMatrix, dof[i], dof[i], rowSum_value[i], INSERT_VALUES); CHKERRV(ierr);      
      }
    }
    
    //For finite element basis greater than 2, to be implemented
    //PetscScalar diag;
    //ierr =MatGetTrace(massMatrix,diag);
    }
    
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
  if(!this->invLumMassMatrixSet())
    this->setInvLumMassMatrix();
  
  const int D = BasisOnMeshType::dim();
  LOG(TRACE)<<"preComputeSystemMatrix" << D << "D";
  
  Mat &systemMatrix = this->data_.systemMatrix_;
  Mat &invLumMassMatrix = this->data_.inversedLumpedMassMatrix_;
  Mat &stiffnessMatrix = this->data_.stiffnessMatrix;
  
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
  
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
solveLinearSystem(Vec &input, Vec &output)
{
  //Au^(t+1)=u^(t)
  
  Mat &systemMatrix=this.systemMatrix;
  
  PetscErrorCode ierr;
  
  // create linear solver context
  std::shared_ptr<Solver::Linear> linearSolver = this->context_.solverManager()->template solver<Solver::Linear>(this->specificSettings_);
  std::shared_ptr<KSP> ksp = linearSolver->ksp();
  
  // set matrix used for linear system and preconditioner to ksp context
  ierr = KSPSetOperators (*ksp, systemMatrix, systemMatrix); CHKERRV(ierr);
  
   // non zero initial values
  PetscScalar scalar = .5;
  ierr = VecSet(output, scalar); CHKERRV(ierr);
  ierr = KSPSetInitialGuessNonzero(*ksp, PETSC_TRUE); CHKERRV(ierr);
  
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