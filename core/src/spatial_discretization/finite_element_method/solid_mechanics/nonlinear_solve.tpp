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
PetscErrorCode nonlinearFunction(SNES snes, Vec u, Vec f, void *context)
{
  T* object = static_cast<T*>(context);
  // note: the input vector u is the same as object->displacement, therefore it does not need to be passed to a method
  
  LOG(DEBUG) << "   (pointer u input: " << u << ")";
  LOG(DEBUG) << "   (pointer displac: " << object->displacements() << ")";
  //assert(u == object->displacements());
  
  if (u != object->displacements())
  {
    object->setDisplacements(u);
    //object->applyDirichletBoundaryConditionsInDisplacements();
  }
  
  // compute the lhs which is the virtual work
  object->computeInternalVirtualWork(f);
  
  // set rows in f to rhs for which dirichlet BC in u is given
  object->applyDirichletBoundaryConditionsInNonlinearFunction(f);
  
  
  
  VLOG(1) << "-- displ in u: " << PetscUtility::getStringVector(u);
  VLOG(1) << "-- computed f: " << PetscUtility::getStringVector(f);
  VLOG(1) << "-- rhs      b: " << PetscUtility::getStringVector(object->rightHandSide());
  
  return 0;
};

/**
 * Evaluates Jacobian matrix
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
PetscErrorCode jacobianFunctionAnalytic(SNES snes, Vec x, Mat jac, Mat b, void *context)
{
  T* object = static_cast<T*>(context);
 
  // note: the input vector x is the same as object->displacement, therefore it does not need to be passed to a method
  
  // compute the tangent stiffness matrix
  object->setStiffnessMatrix();
  
  // store the tangent stiffness matrix in output Vecs
  Mat &tangentStiffnessMatrix = object->tangentStiffnessMatrix();
  //jac = tangentStiffnessMatrix;
  
  assert (tangentStiffnessMatrix == jac);
  assert (tangentStiffnessMatrix == b);
  assert (b == jac);
  
  //MatCopy(tangentStiffnessMatrix,jac,SAME_NONZERO_PATTERN);
  //MatCopy(tangentStiffnessMatrix,b,SAME_NONZERO_PATTERN);
  
  VLOG(1) << "-- computed tangent stiffness matrix analytically: " << PetscUtility::getStringMatrix(jac);
  
  return 0;
};
 

template<typename T>
PetscErrorCode jacobianFunctionFiniteDifferences(SNES snes, Vec x, Mat jac, Mat b, void *context)
{
  T* object = static_cast<T*>(context);
 
  SNESComputeJacobianDefault(snes, x, jac, b, context);
  
  // zero rows and columns for which Dirichlet BC is set 
  object->applyDirichletBoundaryConditionsInStiffnessMatrix(b);
  
  VLOG(1) << "-- computed tangent stiffness matrix by finite differences: " << PetscUtility::getStringMatrix(jac);
  
  return 0;
};

/**
 * Monitor convergence of nonlinear solver
 * 
 *  Input Parameters:
 *   snes  - the SNES context
 *   its   - iteration number
 *   norm  - 2-norm function value (may be estimated)
 *   mctx  - [optional] monitoring context
 */
template<typename T>
PetscErrorCode monitorFunction(SNES snes, PetscInt its, PetscReal norm, void *mctx)
{
  //T* object = static_cast<T*>(mctx);
  VLOG(1) << "  Nonlinear solver: iteration " << its << ", residual norm " << norm;
  return 0;
};
 
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void FiniteElementMethodStiffnessMatrix<
  BasisOnMeshType,
  QuadratureType,
  Term,
  typename BasisOnMeshType::Mesh,
  Equation::isIncompressible<Term>,
  BasisFunction::isNotMixed<typename BasisOnMeshType::BasisFunction>
>::
debug()
{
  LOG(DEBUG) << "------- debug method ------";
 
  //bool useAnalyticJacobian = PythonUtility::getOptionBool(this->specificSettings_, "analyticJacobian", true);
  PetscErrorCode ierr;
  
  //Mat &tangentStiffnessMatrix = this->data_.tangentStiffnessMatrix();
  //Vec &residual = this->data_.residual().values();
  Vec &displacements = this->data_.displacements().values();
  Vec &externalVirtualEnergy = this->data_.externalVirtualEnergy().values();
  
  Vec internalVirtualEnergy;
  VecDuplicate(externalVirtualEnergy, &internalVirtualEnergy);
  
  // zero initial values
  ierr = VecSet(displacements, 0.0); CHKERRV(ierr);
  
  double lx = 1.0;
  double lz = 1.0;
  double tmax = 0.5;
  double that = 0.5;
  
  const double c0 = SEMT::Parameter<0>::get_value();
  const double c1 = SEMT::Parameter<1>::get_value();
  assert(c1 == 0);
  
  //double lambdaValue = pow(1 + tmax/(2*c0), 3);
  
  // constant load per surface area, t = that, but total load varying
  double lambdaValue = 0.30285343213869*pow(c0,(-0.5))*
  pow(18.0*pow(c0,1.5) + 2.44948974278318*pow((54.0*pow(c0,3) 
  - pow(that,3)),0.5),-0.333333333333333)
  *(1.81712059283214*that + pow((18.0*pow(c0,1.5) 
  + 2.44948974278318*pow((54.0*pow(c0,3) - pow(that,3)),0.5)),0.666666666666667));
  
  LOG(DEBUG) << "c0: " << c0 << ", c1: " << c1 << ", tmax: " << tmax << ", lambdaValue: " << lambdaValue;

  LOG(DEBUG) << "expected F: diag( " << lambdaValue << ", " << 1./sqrt(lambdaValue) << ", " << 1./sqrt(lambdaValue) << ")";
  LOG(DEBUG) << "expected C: diag( " << lambdaValue*lambdaValue << ", " << 1./lambdaValue << ", " << 1./lambdaValue << ")";
  LOG(DEBUG) << "expected I1: " << lambdaValue*lambdaValue + 2/lambdaValue;
  LOG(DEBUG) << "expected S: diag( " << 2*c0 + 4*c1/lambdaValue << ", " << -2*c1/lambdaValue << ", " << -2*c1/lambdaValue <<  " )";
  
  double s11 = 2*c0 + 4*c1/lambdaValue;
  
  LOG(DEBUG) << "expected W[9]: " << 1./4*lz*lz*lambdaValue*s11 << " = " << 1./4*tmax;
  
  std::vector<double> knownValues = {
   0.0, 0.0, 0.0, 
   (lambdaValue-1.0)*lx, 0.0, 0.0,
   0.0, (1./sqrt(lambdaValue) - 1.0)*lz, 0.0, 
   (lambdaValue-1.0)*lx, (1./sqrt(lambdaValue) - 1.0)*lz, 0.0,
   0.0, 0.0, (1./sqrt(lambdaValue) - 1.0)*lz, 
   (lambdaValue-1.0)*lx, 0.0, (1./sqrt(lambdaValue) - 1.0)*lz,
   0.0, (1./sqrt(lambdaValue) - 1.0)*lz, (1./sqrt(lambdaValue) - 1.0)*lz, 
   (lambdaValue-1.0)*lx, (1./sqrt(lambdaValue) - 1.0)*lz, (1./sqrt(lambdaValue) - 1.0)*lz
  };
  PetscUtility::setVector(knownValues, displacements);
  
  // set prescribed Dirchlet BC displacements values
  //applyDirichletBoundaryConditionsInDisplacements();
  
  VLOG(1) << "displacements values: " << PetscUtility::getStringVector(displacements);
  VLOG(1) << "W_ext: " << PetscUtility::getStringVector(externalVirtualEnergy);
  
  this->computeInternalVirtualWork(internalVirtualEnergy);
  VLOG(1) << "W_int: " << PetscUtility::getStringVector(internalVirtualEnergy);
  
  std::vector<double> wExt, wInt;
  PetscUtility::getVectorEntries(externalVirtualEnergy, wExt);
  PetscUtility::getVectorEntries(internalVirtualEnergy, wInt);
  
  std::stringstream s;
  s << std::endl << "no.  Wext  Wint" << std::endl;
  for (int i=0; i<wExt.size(); i++)
  {
    s << i << "   " << wExt[i] << "   " << wInt[i] << std::endl;
  }
  LOG(DEBUG) << s.str();
  
  LOG(DEBUG) << "------- debug method end, exit------";
  exit(0);
  
}

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
  
  debug();
  exit(0);
  
  bool useAnalyticJacobian = PythonUtility::getOptionBool(this->specificSettings_, "analyticJacobian", true);
  PetscErrorCode ierr;
  
  Mat &tangentStiffnessMatrix = this->data_.tangentStiffnessMatrix();
  Vec &residual = this->data_.residual().values();
  Vec &displacements = this->data_.displacements().values();
  Vec &externalVirtualEnergy = this->data_.externalVirtualEnergy().values();
  
  LOG(DEBUG) << "residual: " << residual;
  LOG(DEBUG) << "displacements: " << displacements;
  
  // create nonlinear solver PETSc context (snes)
  std::shared_ptr<Solver::Nonlinear> nonlinearSolver = this->context_.solverManager()->template solver<Solver::Nonlinear>(this->specificSettings_);
  std::shared_ptr<SNES> snes = nonlinearSolver->snes();
  std::shared_ptr<KSP> ksp = nonlinearSolver->ksp();
  
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
  PetscErrorCode (*callbackJacobianAnalytic)(SNES, Vec, Mat, Mat, void *) = *jacobianFunctionAnalytic<ThisClass>;
  PetscErrorCode (*callbackJacobianFiniteDifferences)(SNES, Vec, Mat, Mat, void *) = *jacobianFunctionFiniteDifferences<ThisClass>;
  PetscErrorCode (*callbackMonitorFunction)(SNES, PetscInt, PetscReal, void *) = *monitorFunction<ThisClass>;
  
  // set function
  ierr = SNESSetFunction(*snes, residual, callbackNonlinearFunction, this); CHKERRV(ierr);
  
  // set jacobian 
  if (useAnalyticJacobian)
  {
    ierr = SNESSetJacobian(*snes, tangentStiffnessMatrix, tangentStiffnessMatrix, callbackJacobianAnalytic, this); CHKERRV(ierr);
    LOG(DEBUG) << "Use analytical jacobian";
  }
  else 
  {
    // set function to compute jacobian from finite differences
    ierr = SNESSetJacobian(*snes, tangentStiffnessMatrix, tangentStiffnessMatrix, callbackJacobianFiniteDifferences, this); CHKERRV(ierr);
    LOG(DEBUG) << "Use Finite-Differences approximation for jacobian";
  }

  // set monitor function  
  ierr = SNESMonitorSet(*snes, callbackMonitorFunction, this, NULL); CHKERRV(ierr);

  // zero initial values
  ierr = VecSet(displacements, 0.0); CHKERRV(ierr);
  
  // set prescribed Dirchlet BC displacements values
  applyDirichletBoundaryConditionsInDisplacements();
  
  VLOG(1) << "Dirichlet BC: indices: " << this->dirichletIndices_ << ", values: " << dirichletValues_ << ", rhsValues: " << rhsValues_;
  VLOG(1) << "initial values: " << PetscUtility::getStringVector(displacements);

  if (!useAnalyticJacobian)
  {
     LOG(DEBUG) << "compute analytical jacobian to initialize non-zero structure of matrix";
     callbackJacobianAnalytic(*snes, displacements, tangentStiffnessMatrix, tangentStiffnessMatrix, this);
  }
  
  LOG(DEBUG) << "------------------  start solve  ------------------";
  
  // solve the system nonlinearFunction(displacements) = externalVirtualEnergy
  // not sure if externalVirtualEnergy and displacements have to be different vectors from the ones used in the provided functions
  ierr = SNESSolve(*snes, externalVirtualEnergy, displacements); CHKERRV(ierr);
  
  int numberOfIterations = 0;
  PetscReal residualNorm = 0.0;
  ierr = SNESGetIterationNumber(*snes, &numberOfIterations); CHKERRV(ierr);
  ierr = SNESGetFunctionNorm(*snes, &residualNorm); CHKERRV(ierr);
  
  SNESConvergedReason convergedReason;
  KSPConvergedReason kspConvergedReason;
  ierr = SNESGetConvergedReason(*snes, &convergedReason); CHKERRV(ierr);
  ierr = KSPGetConvergedReason(*ksp, &kspConvergedReason); CHKERRV(ierr);
  
  LOG(INFO) << "Solution done in " << numberOfIterations << " iterations, residual norm " << residualNorm 
    << ": " << PetscUtility::getStringNonlinearConvergedReason(convergedReason) << ", " 
    << PetscUtility::getStringLinearConvergedReason(kspConvergedReason);
  
  // update geometry field from displacements, w = alpha*x+y
  VecWAXPY(this->data_.geometryActual().values(), 1.0, this->data_.geometryReference().values(), displacements);
}
  
};