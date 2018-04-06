#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_incompressible.h"

#include <Python.h>  // this has to be the first included header 
#include <iostream>
#include <petscsnes.h>
#include <petscksp.h>

#include "easylogging++.h"
#include "solver/nonlinear.h"
#include "data_management/finite_elements_solid_mechanics.h"

namespace SpatialDiscretization
{
 
/** 
 * Nonlinear function F that gets solved by PETSc, solve for x such that F(x) = b
 *  Input Parameters:
 *  snes - the SNES context
 *  x    - input vector
 *  context  - optional user-defined context
 *
 *  Output Parameter:
 *  f - function vector
 */
template<typename T>
PetscErrorCode nonlinearFunction(SNES snes, Vec x, Vec f, void *context)
{
  T* object = static_cast<T*>(context);
 
  // save the current displacement values
  object->setDisplacements(x);
  
  // compute the lhs which is the virtual work
  object->computeInternalVirtualWork(f);
  
  // set rows in f to rhs for which dirichlet BC in u is given
  object->applyDirichletBoundaryConditionsInNonlinearFunction(f);
  
  return 0;
};

/**
 * FormJacobian1 - Evaluates Jacobian matrix
 *  Input Parameters:
 *  snes - the SNES context
 *  x    - input vector
 *  context  - optional user-defined context
 *
 *  Output Parameters:
 *  jac - Jacobian matrix
 *  b   - optionally different preconditioning matrix
 */
template<typename T>
PetscErrorCode jacobianFunction(SNES snes, Vec x, Mat jac, Mat b, void *context)
{
  T* object = static_cast<T*>(context);
 
  // save the current displacement values
  object->setDisplacements(x);
  
  // compute the tangent stiffness matrix
  object->setStiffnessMatrix();
  
  // store the tangent stiffness matrix in output Vecs
  Mat &tangentStiffnessMatrix = object->tangentStiffnessMatrix();
  jac = tangentStiffnessMatrix;
  b = jac;
  return 0;
};
 
// general implementation of nonlinear solving for solid mechanics (here for penalty formulation)
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
solve()
{
  LOG(TRACE) << "FiniteElementMethod::solve (nonlinear)";
  
  PetscErrorCode ierr;
  
  Mat &tangentStiffnessMatrix = this->data_.tangentStiffnessMatrix();
  Vec &residual = this->data_.residual().values();
  Vec &displacements = this->data_.displacements().values();
  Vec &externalVirtualEnergy = this->data_.externalVirtualEnergy().values();
  
  // create nonlinear solver PETSc context (snes)
  std::shared_ptr<Solver::Nonlinear> nonlinearSolver = this->context_.solverManager()->template solver<Solver::Nonlinear>(this->specificSettings_);
  std::shared_ptr<SNES> snes = nonlinearSolver->snes();
  
  assert(snes != nullptr);
  
  typedef FiniteElementMethodStiffnessMatrix<
    BasisOnMeshType,
    QuadratureType,
    Term,
    typename BasisOnMeshType::Mesh,
    Equation::isIncompressible<Term>,
    BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
  > ThisClass;
  
  PetscErrorCode (*callbackNonlinearFunction)(SNES, Vec, Vec, void *) = *nonlinearFunction<ThisClass>;
  PetscErrorCode (*callbackJacobian)(SNES, Vec, Mat, Mat, void *) = *jacobianFunction<ThisClass>;
  
  // set function
  SNESSetFunction(*snes, residual, callbackNonlinearFunction, this);
  
  // set jacobian 
  SNESSetJacobian(*snes, tangentStiffnessMatrix, tangentStiffnessMatrix, callbackJacobian, this);
  
  // zero initial values
  ierr = VecSet(displacements, 0.0); CHKERRV(ierr);
  
  // solve the system nonlinearFunction(displacements) = externalVirtualEnergy
  // not sure if externalVirtualEnergy and displacements have to be different vectors from the ones used in the provided functions
  ierr = SNESSolve(*snes, externalVirtualEnergy, displacements); CHKERRV(ierr);
  
  int numberOfIterations = 0;
  PetscReal residualNorm = 0.0;
  ierr = SNESGetIterationNumber(*snes, &numberOfIterations); CHKERRV(ierr);
  ierr = SNESGetFunctionNorm(*snes, &residualNorm); CHKERRV(ierr);
  
  SNESConvergedReason convergedReason;
  ierr = SNESGetConvergedReason(*snes, &convergedReason); CHKERRV(ierr);
  
  LOG(INFO) << "Solution done in " << numberOfIterations << " iterations, residual norm " << residualNorm 
    << ": " << PetscUtility::getStringNonlinearConvergedReason(convergedReason);
  
  // update geometry field from displacements, w = alpha*x+y
  VecWAXPY(this->data_.geometryActual().values(), 1.0, this->data_.geometryReference().values(), displacements);
}
  
};