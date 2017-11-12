#include "spatial_discretization/finite_element_method.h"

#include <iostream>
#include <petscksp.h>
#include <vector>
#include <numeric>

#include <Python.h>
#include "easylogging++.h"

#include "control/types.h"
#include "control/python_utility.h"
#include <control/petsc_utility.h>


namespace SpatialDiscretization
{
  
template<class MeshType, class BasisFunctionType>
FiniteElementMethodBaseTimeStepping<MeshType, BasisFunctionType>::
FiniteElementMethodBaseTimeStepping(DihuContext &context)
  : FiniteElementMethodBase<MeshType, BasisFunctionType>(context)
{
}

template<typename MeshType, typename BasisFunctionType, typename Term>
FiniteElementMethod<MeshType, BasisFunctionType, Term, Equation::hasLaplaceOperatorWithTimeStepping<Term>>::
FiniteElementMethod(DihuContext &context) 
  : FiniteElementMethodBaseTimeStepping<MeshType, BasisFunctionType>(context)
{
}

template<typename Mesh, typename BasisFunction>
void FiniteElementMethodBaseTimeStepping<Mesh, BasisFunction>::
initialize()
{
  this->setStiffnessMatrix();
  this->createRhsDiscretizationMatrix();
  this->data_.finalAssembly();
}

template<class MeshType, class BasisFunctionType>
void FiniteElementMethodBaseTimeStepping<MeshType, BasisFunctionType>::
createRhsDiscretizationMatrix()
{
}

template<class MeshType, class BasisFunctionType>
void FiniteElementMethodBaseTimeStepping<MeshType, BasisFunctionType>::
recoverRightHandSide(Vec &result)
{
  // dmatrix * f_strong = rhs_weak
  Vec &rhs = this->data_.rightHandSide();   // rhs in weak formulation
  Mat &dmatrix = this->data_.discretizationMatrix();
  
  LOG(DEBUG) << "recoverRightHandSide";
  
  PetscErrorCode ierr;
  
  // create linear solver context
  KSP ksp; 
  ierr = KSPCreate (PETSC_COMM_WORLD, &ksp); CHKERRV(ierr);  
  
  // set matrix used for linear system and preconditioner to ksp context
  ierr = KSPSetOperators (ksp, dmatrix, dmatrix); CHKERRV(ierr);
  
  // extract preconditioner context
  PC pc;
  ierr = KSPGetPC (ksp, &pc); CHKERRV(ierr);
  
  // set preconditioner type
  ierr = PCSetType (pc, PCJACOBI); CHKERRV(ierr);
  
  // set solution tolerances
  double relativeTolerance = PythonUtility::getOptionDouble(this->specificSettings_, "relativeTolerance", 1e-5, PythonUtility::Positive);

  //                            relative tol,      absolute tol,  diverg tol.,   max_iterations
  ierr = KSPSetTolerances (ksp, relativeTolerance, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRV(ierr);

  // non zero initial values
  PetscScalar scalar = .5;
  ierr = VecSet(result, scalar); CHKERRV(ierr);
  ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERRV(ierr);
  
  // solve the system
  ierr = KSPSolve(ksp, rhs, result); CHKERRV(ierr);
  
  int numberOfIterations = 0;
  PetscReal residualNorm = 0.0;
  ierr = KSPGetIterationNumber(ksp, &numberOfIterations); CHKERRV(ierr);
  ierr = KSPGetResidualNorm(ksp, &residualNorm); CHKERRV(ierr);
  
  LOG(INFO) << "Rhs recovered in " << numberOfIterations << " iterations, residual norm " << residualNorm;
}

template<class MeshType, class BasisFunctionType>
void FiniteElementMethodBaseTimeStepping<MeshType, BasisFunctionType>::
evaluateTimesteppingRightHandSide(Vec &input, Vec &output)
{
  Mat &stiffnessMatrix = this->data_.stiffnessMatrix();
  Vec &rhs = this->data_.rightHandSide();
  
  // compute rhs = stiffnessMatrix*input
  MatMult(stiffnessMatrix, input, rhs);
  
  recoverRightHandSide(output);
  
  this->data_.print();  
}

} // namespace SpatialDiscretization