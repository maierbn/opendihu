#include "model_order_reduction/pod.h"

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

namespace ModelOrderReduction
{


template<typename DiscretizableInTimeType>
PODBase<DiscretizableInTimeType>::
PODBase(DihuContext context) :
  DiscretizableInTime(SolutionVectorMapping(true)), problem_(context), context_(context)
{
}

template<typename DiscretizableInTimeType>
void PODBase<DiscretizableInTimeType>::
initialize()
{
  // es gibt 2 versch. initialize methoden: mit und ohne arg.
  problem_.initialize();
}

template<typename DiscretizableInTimeType>
void PODBase<DiscretizableInTimeType>::
preComputeSystemMatrix(double timeStepWidth)
{
  // redSysMatr = V^t * SysMatr * V MatMatMatMult
  problem_.preComputeSystemMatrix(timeStepWidth);
  problem_->data.systemMatrix();
  MatCreateTranspose(V, &Vt);
  // data protected, aber data() Methode public
  Mat &systemMatrix = problem_->data.systemMatrix();
  PetscErrorCode ierr;
  ierr = MatMatMatMult(Vt, systemMatrix, V, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &redSysMat);
}

template<typename DiscretizableInTimeType>
void PODBase<DiscretizableInTimeType>::
solveLinearSystem(Vec &input, Vec &output)
{
  // V^t * input = z(t)
  // solve imlicit z(t+1) = redSysMatr * z(t)
  PetscErrorCode ierr;
  Vec temp;
  ierr = MatMult(Vt, input, temp);
  
  // create linear solver context
  // problem: wie wird das hier richtig aufgerufen
  // context private (selber erstellen?)
  // problem: schon erstellten solver aufrufen, werte in settings korrekt für diesen Fall?
  std::shared_ptr<Solver::Linear> linearSolver = problem_->context_.solverManager()->template solver<Solver::Linear>(problem_->specificSettings_);
  std::shared_ptr<KSP> ksp = linearSolver->ksp();
  
  // set matrix used for linear system and preconditioner to ksp context
  ierr = KSPSetOperators (*ksp, redSysMat, redSysMat); CHKERRV(ierr);
  
  //use the default initial guess (zero) by the PETSC because the input and output are the same by the implicit method   
  // solve the system
  ierr = KSPSolve(*ksp, temp, output); CHKERRV(ierr);
  
  int numberOfIterations = 0;
  PetscReal residualNorm = 0.0;
  ierr = KSPGetIterationNumber(*ksp, &numberOfIterations); CHKERRV(ierr);
  ierr = KSPGetResidualNorm(*ksp, &residualNorm); CHKERRV(ierr);
  
}

//! timestepping rhs function f of equation u_t = f(u,t)
template<typename DiscretizableInTimeType>
void PODBase<DiscretizableInTimeType>::
evaluateTimesteppingRightHandSide(Vec &input, Vec &output, int timeStepNo, double currentTime)
{
  // call the method of the underlying problem
  problem_.evaluateTimesteppingRightHandSide(input, output, timeStepNo, currentTime);
}

template<typename DiscretizableInTimeType>
void POD<DiscretizableInTimeType, LinearPart>::
evaluateTimesteppingRightHandSide(Vec &input, Vec &output, int timeStepNo, double currentTime)
{
  PODBase<DiscretizableInTimeType>::evaluateTimesteppingRightHandSide(input, output, timeStepNo, currentTime);
};

//! get the number of degrees of freedom per node which is 1 by default
template<typename DiscretizableInTimeType>
int PODBase<DiscretizableInTimeType>::
nComponentsNode()
{
  return problem_.nComponentsNode();
}

//! set initial values and return true or don't do anything and return false
template<typename DiscretizableInTimeType>
bool PODBase<DiscretizableInTimeType>::
setInitialValues(Vec &initialValues)
{
  problem_.setInitialValues(initialValues);
}

//! return whether the object has a specified mesh type and is not independent of the mesh type
template<typename DiscretizableInTimeType>
bool PODBase<DiscretizableInTimeType>::
knowsMeshType()
{
  return problem_.knowsMeshType();
}
template<typename DiscretizableInTimeType>
std::shared_ptr<Mesh::Mesh> PODBase<DiscretizableInTimeType>::
mesh()
{
  return problem_.mesh();
}

};
