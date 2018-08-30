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
evaluateTimesteppingRightHandSideImplicit(Vec &input, Vec &output, int timeStepNo, double currentTime)
{
  LOG(TRACE) << "evaluateTimesteppingRightHandSideImplicit";
  
  // this method computes output = systemMatrix^{-1}*input with systemMatrix = (I - dt*M^{-1}K)

  // solve the linear system
  solveLinearSystem(input, output);
  
  VLOG(1) << PetscUtility::getStringVector(output);
}
/*
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
precomputeSystemMatrix1()
{
  const int D = BasisOnMeshType::dim();
  LOG(TRACE) << "precomputeSystemMatrix1 " << D << "D";
  
  //modify the stiffness matrix to use it as the system matrix
  this->data_.setSystemMatrix(this->data_.stiffnessMatrix());
  Mat &massMatrix = this->data_.massMatrix().valuesGlobal();
  

  //Scale is 1.0 because the mass matrix is scaled with the time step by the initialization
  PetscScalar scale = 1.0;
  // compute systemMatrix += scale*massMatrix
  PetscErrorCode ierr;
  ierr = MatAXPY(this->data_.systemMatrix()->valuesGlobal(), scale, massMatrix, SAME_NONZERO_PATTERN); CHKERRV(ierr);
  this->data_.systemMatrix()->assembly(MAT_FINAL_ASSEMBLY);
}
*/
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
setInverseLumpedMassMatrix()
{
  LOG(TRACE) << "setInverseLumpedMassMatrix";

  std::shared_ptr<BasisOnMeshType> mesh = std::static_pointer_cast<BasisOnMeshType>(this->data_.mesh());
  Mat &inverseLumpedMassMatrix = this->data_.inverseLumpedMassMatrix()->globalValues();
  Mat &massMatrix = this->data_.massMatrix()->globalValues();
  
  PetscErrorCode ierr;
  
  PetscInt nRows, nColumns;
  ierr = MatGetSize(massMatrix,&nRows,&nColumns); CHKERRV(ierr);
  VLOG(1) << "massMatrix nRows " << nRows << " nColumns " << nColumns;
     
  std::shared_ptr<PartitionedPetscVec<BasisOnMeshType,1>> rowSum = std::make_shared<PartitionedPetscVec<BasisOnMeshType,1>>(mesh->meshPartition(), "rowSum");

  // In case of linear and bilinear basis functions
  // store the sum of each row of the matrix in the vector rowSum
  ierr = MatGetRowSum(massMatrix, rowSum->valuesGlobal()); CHKERRV(ierr);

  // for the inverse matrix, replace each entry in rowSum by its reciprocal
  ierr = VecReciprocal(rowSum); CHKERRV(ierr);

  // set the values on the diagonal
  ierr = MatDiagonalSet(inverseLumpedMassMatrix, rowSum->valuesGlobal(), INSERT_VALUES); CHKERRV(ierr);

  this->data_.inverseLumpedMassMatrix()->assembly(MAT_FINAL_ASSEMBLY);

  VLOG(2) << *this->data_.inverseLumpedMassMatrix();
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
setSystemMatrix(double timeStepWidth)
{
  LOG(TRACE) << "setSystemMatrix(timeStepWidth=" << timeStepWidth << ")";

  // compute the system matrix (I - dt*M^{-1}K) where M^{-1} is the lumped mass matrix

  Mat &inverseLumpedMassMatrix = this->data_.inverseLumpedMassMatrix()->globalValues();
  Mat &stiffnessMatrix = this->data_.stiffnessMatrix()->globalValues();
  Mat systemMatrix;
  
  PetscErrorCode ierr;
       
  // compute systemMatrix = M^{-1}K
  // the result matrix is created by MatMatMult
  ierr = MatMatMult(inverseLumpedMassMatrix, stiffnessMatrix, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &systemMatrix);
  this->data_.initializeSystemMatrix(systemMatrix);
  
  // scale systemMatrix by -dt, systemMatrix = -dt*M^{-1}K
  ierr = MatScale(this->data_.systemMatrix().valuesGlobal(), -timeStepWidth); CHKERRV(ierr);

  // add 1 on the diagonal: systemMatrix = I - dt*M^{-1}K
  ierr = MatShift(this->data_.systemMatrix().valuesGlobal(), 1.0); CHKERRV(ierr);

  this->data_.systemMatrix().assembly(MAT_FINAL_ASSEMBLY);

  VLOG(1) << *this->data_.systemMatrix();
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
solveLinearSystem(Vec &input, Vec &output)
{
  // solve systemMatrix*output = input for output
  Mat &systemMatrix = this->data_.systemMatrix()->valuesGlobal();

  PetscErrorCode ierr;
  PetscUtility::checkDimensionsMatrixVector(systemMatrix, input);

  // set matrix used for linear system and preconditioner to ksp context
  assert(this->ksp_);
  ierr = KSPSetOperators(*ksp_, systemMatrix, systemMatrix); CHKERRV(ierr);

  // solve the system
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

} // namespace SpatialDiscretization
