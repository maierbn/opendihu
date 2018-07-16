#include "spatial_discretization/finite_element_method/05_timestepping.h"

#include <Python.h>
#include <iostream>
#include <sstream>
#include <petscksp.h>
#include <vector>
#include <numeric>
#include <omp.h>

#include "easylogging++.h"

#include "control/types.h"
#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "solver/solver_manager.h"
#include "solver/linear.h"


namespace SpatialDiscretization
{

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
FiniteElementMethodTimeStepping(DihuContext context)
  : FiniteElementMethodBaseRhs<BasisOnMeshType, QuadratureType, Term>(context),
  DiscretizableInTime(SolutionVectorMapping(true)), linearSolver_(nullptr), ksp_(nullptr)
{
  // the solutionVectorMapping_ object stores the information which range of values of the solution will be further used
  // in methods that use the result of this method, e.g. in operator splittings. Since there are no internal values
  // in this FEM, set the range to all values.
  solutionVectorMapping_.setOutputRange(0, this->data_.mesh()->nLocalNodes());
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
initialize()
{
  this->data_.initialize();
  this->setStiffnessMatrix();
  this->setMassMatrix();
  this->data_.finalAssembly();
  initializeLinearSolver();

  if (this->outputWriterManager_.hasOutputWriters())
  {
    LOG(WARNING) << "You have specified output writers for a FiniteElementMethod which is used for a time stepping problem. "
      "The output will not contain any solution data, only geometry. Probably you want to get output from the time stepping scheme, then define the output writers there.";
  }
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
reset()
{
  linearSolver_ = nullptr;
  ksp_ = nullptr;
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
setRankSubset(Partition::RankSubset rankSubset)
{
  FiniteElementMethodBase<BasisOnMeshType, QuadratureType, Term>::setRankSubset(rankSubset);
}
  
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
initializeLinearSolver()
{ 
  if (linearSolver_ == nullptr)
  {
    std::stringstream s;
    s << "[" << omp_get_thread_num() << "/" << omp_get_num_threads() << "]";
    LOG(DEBUG) << s.str() << ": FiniteElementMethodTimeStepping: initialize linearSolver";
    
    // retrieve linear solver
    linearSolver_ = this->context_.solverManager()->template solver<Solver::Linear>(this->specificSettings_);
    ksp_ = linearSolver_->ksp();
  }
  else 
  {
    std::stringstream s;
    s << "[" << omp_get_thread_num() << "/" << omp_get_num_threads() << "]";
    VLOG(2) << s.str() << ": linearSolver_ already set";
    
  }
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
constexpr int FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
nComponents()
{
  return 1;   // this may be different for structural mechanics
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
recoverRightHandSideStrongForm(Vec &result)
{
  // massMatrix * f_strong = rhs_weak
  Vec &rhs = this->data_.rightHandSide().values();   // rhs in weak formulation
  std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> massMatrix = this->data_.massMatrix();

  PetscErrorCode ierr;

  // create linear solver context
  VLOG(1) << omp_get_thread_num() << ": recoverRightHandSideStrongForm";
  initializeLinearSolver();
  
  // set matrix used for linear system and preconditioner to ksp context
  ierr = KSPSetOperators (*ksp_, massMatrix, massMatrix); CHKERRV(ierr);

  // non zero initial values
  PetscScalar scalar = .5;
  ierr = VecSet(result, scalar); CHKERRV(ierr);
  ierr = KSPSetInitialGuessNonzero(*ksp_, PETSC_TRUE); CHKERRV(ierr);

  // solve the system
  ierr = KSPSolve(*ksp_, rhs, result); CHKERRV(ierr);

  int numberOfIterations = 0;
  PetscReal residualNorm = 0.0;
  ierr = KSPGetIterationNumber(*ksp_, &numberOfIterations); CHKERRV(ierr);
  ierr = KSPGetResidualNorm(*ksp_, &residualNorm); CHKERRV(ierr);

  //LOG(INFO) << "Rhs recovered in " << numberOfIterations << " iterations, residual norm " << residualNorm;
  VLOG(1) << "Rhs recovered in " << numberOfIterations << " iterations, residual norm " << residualNorm;
}


template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
checkDimensions(std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> stiffnessMatrix, Vec &input)
{
#ifndef NDEBUG
  int nRows, nColumns;
  MatGetSize(stiffnessMatrix, &nRows, &nColumns);
  int nEntries;
  VecGetSize(input, &nEntries);
  if (nColumns != nEntries)
  {
    LOG(ERROR) << "Stiffness matrix dimension " << nRows << "x" << nColumns << " does not match input vector (" << nEntries << ")!";
  }
  assert(nColumns == nEntries);
#endif
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
evaluateTimesteppingRightHandSide(Vec &input, Vec &output, int timeStepNo, double currentTime)
{
  std::shared_ptr<PartitionedPetscMat<BasisOnMeshType>> stiffnessMatrix = this->data_.stiffnessMatrix();
  Vec &rhs = this->data_.rightHandSide().values();

  // check if matrix and vector sizes match
  checkDimensions(stiffnessMatrix, input);

  // compute rhs = stiffnessMatrix*input
  MatMult(stiffnessMatrix, input, rhs);
 
  // compute output = massMatrix^{-1}*rhs
  recoverRightHandSideStrongForm(output);

  this->data_.print();
  this->outputWriterManager_.writeOutput(this->data_, timeStepNo, currentTime);
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
bool FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
knowsMeshType()
{
  return true;
}

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
std::shared_ptr<Mesh::Mesh> FiniteElementMethodTimeStepping<BasisOnMeshType, QuadratureType, Term>::
mesh()
{
  return FiniteElementMethodBase<BasisOnMeshType, QuadratureType, Term>::mesh();
}

} // namespace SpatialDiscretization