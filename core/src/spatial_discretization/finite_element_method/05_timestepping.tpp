#include "spatial_discretization/finite_element_method/05_timestepping.h"

#include <Python.h>
#include <iostream>
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
  
template<typename BasisOnMeshType, typename IntegratorType, typename Term>
FiniteElementMethodTimeStepping<BasisOnMeshType, IntegratorType, Term>::
FiniteElementMethodTimeStepping(const DihuContext &context)
  : FiniteElementMethodBaseRhs<BasisOnMeshType, IntegratorType, Term>(context),
  DiscretizableInTime(SolutionVectorMapping(true))
{
  // the solutionVectorMapping_ object stores the information which range of values of the solution will be further used 
  // in methods that use the result of this method, e.g. in operator splittings. Since there are no internal values
  // in this FEM, set the range to all values.
  solutionVectorMapping_.setOutputRange(0, this->data_.mesh()->nNodes());
}

template<typename BasisOnMeshType, typename IntegratorType, typename Term>
void FiniteElementMethodTimeStepping<BasisOnMeshType, IntegratorType, Term>::
initialize()
{
  this->data_.initialize();
  this->setStiffnessMatrix();
  this->setRhsDiscretizationMatrix();
  this->data_.finalAssembly();
}

template<typename BasisOnMeshType, typename IntegratorType, typename Term>
void FiniteElementMethodTimeStepping<BasisOnMeshType, IntegratorType, Term>::
recoverRightHandSide(Vec &result)
{
  // discretizationMatrix * f_strong = rhs_weak
  Vec &rhs = this->data_.rightHandSide().values();   // rhs in weak formulation
  Mat &discretizationMatrix = this->data_.discretizationMatrix();
  
  PetscErrorCode ierr;
  
  // create linear solver context
  std::shared_ptr<Solver::Linear> linearSolver = this->context_.solverManager()->template solver<Solver::Linear>(this->specificSettings_);
  std::shared_ptr<KSP> ksp = linearSolver->ksp();
  
  // set matrix used for linear system and preconditioner to ksp context
  ierr = KSPSetOperators (*ksp, discretizationMatrix, discretizationMatrix); CHKERRV(ierr);
  
  // non zero initial values
  PetscScalar scalar = .5;
  ierr = VecSet(result, scalar); CHKERRV(ierr);
  ierr = KSPSetInitialGuessNonzero(*ksp, PETSC_TRUE); CHKERRV(ierr);
  
  // solve the system
  ierr = KSPSolve(*ksp, rhs, result); CHKERRV(ierr);
  
  int numberOfIterations = 0;
  PetscReal residualNorm = 0.0;
  ierr = KSPGetIterationNumber(*ksp, &numberOfIterations); CHKERRV(ierr);
  ierr = KSPGetResidualNorm(*ksp, &residualNorm); CHKERRV(ierr);
  
  //LOG(INFO) << "Rhs recovered in " << numberOfIterations << " iterations, residual norm " << residualNorm;
}

template<typename BasisOnMeshType, typename IntegratorType, typename Term>
void FiniteElementMethodTimeStepping<BasisOnMeshType, IntegratorType, Term>::
evaluateTimesteppingRightHandSide(Vec &input, Vec &output, int timeStepNo, double currentTime)
{
  Mat &stiffnessMatrix = this->data_.stiffnessMatrix();
  Vec &rhs = this->data_.rightHandSide().values();
  
  // compute rhs = stiffnessMatrix*input
  MatMult(stiffnessMatrix, input, rhs);
  
  recoverRightHandSide(output);
  
  this->data_.print();
}

template<typename BasisOnMeshType, typename IntegratorType, typename Term>
bool FiniteElementMethodTimeStepping<BasisOnMeshType, IntegratorType, Term>::
knowsMeshType()
{
  return true;
}
  
} // namespace SpatialDiscretization