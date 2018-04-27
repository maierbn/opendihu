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
  
  //LOG(DEBUG) << "   nonlinear function, pointer u input: " << u;
  //LOG(DEBUG) << "   nonlinear function, pointer displac: " << object->displacements();
  //assert(u == object->displacements());
      
  // note: if the input vector u is the same as object->displacement, it does not need to be passed to a method
  if (u != object->displacements())
  {
    object->setDisplacements(u);
  }
  
  object->applyDirichletBoundaryConditionsInDisplacements();
  
  // compute the lhs which is the virtual work
  object->computeInternalMinusExternalVirtualWork(f);
  
  // set rows in f to 0 for which dirichlet BC in u is given
  object->applyDirichletBoundaryConditionsInNonlinearFunction(f);
  
  LOG(DEBUG) << "--            displ input u with BC: " << PetscUtility::getStringVector(u);
  //LOG(DEBUG) << "--            displ input u with BC: " << PetscUtility::getStringVector(object->displacements());
  LOG(DEBUG) << "-- computed f=dW_int-dW_ext with BC: " << PetscUtility::getStringVector(f);
  
  // compute and output function norm
  PetscReal functionNorm;
  VecNorm(f, NORM_2, &functionNorm);
  LOG(DEBUG) << "function norm: " << functionNorm;
  
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
 
  // note: if the input vector x is the same as object->displacement, it does not need to be passed to a method
  if (x != object->displacements())
  {
    object->setDisplacements(x);
  }
  
  // compute the tangent stiffness matrix
  object->computeAnalyticalStiffnessMatrix(jac);
  
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
  LOG(DEBUG) << "  Nonlinear solver: iteration " << its << ", residual norm " << norm;
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
  Vec &externalVirtualWork = this->data_.externalVirtualWork().values();
  
  Vec internalVirtualWork;
  VecDuplicate(externalVirtualWork, &internalVirtualWork);
  
  // zero initial values
  ierr = VecSet(displacements, 0.0); CHKERRV(ierr);
 
  if (BasisOnMeshType::dim() == 3)
  {
    double lx = 1.0;
    double lz = 1.0;
    double tmax = 0.5;
    
    const double c0 = SEMT::Parameter<0>::get_value();
    const double c1 = SEMT::Parameter<1>::get_value();
    assert(c1 == 0);
    
    double lambdaValue = 1.0;
    
    // constant total load
    
    // analytical solution from equilibrium
    lambdaValue = MathUtility::sqr(tmax)/(6*c0*MathUtility::sqr(lz)*pow((108*pow(c0,3)*pow(lz,6) + pow(tmax,3) + sqrt(-pow(tmax,6)
    + pow((108*pow(c0,3)*pow(lz,6) + pow(tmax,3)),2))),1./3.))
    + tmax/(6*c0*MathUtility::sqr(lz))
    + pow(108*pow(c0,3)*pow(lz,6) + pow(tmax,3) + sqrt(-pow(tmax,6) + MathUtility::sqr(108*pow(c0,3)*pow(lz,6) + pow(tmax,3))),1./3.)
    /(6*c0*MathUtility::sqr(lz));
    
    // solution from penalty formulation 
    lambdaValue = MathUtility::sqr(tmax)/(6*c0*MathUtility::sqr(lz)*pow((108*pow(c0,3)*pow(lz,6) + pow(tmax,3)
      + sqrt(-pow(tmax,6) + pow((108*pow(c0,3)*pow(lz,6)
      + pow(tmax,3)),2))),1./3.)) + tmax/(6*c0*MathUtility::sqr(lz)) + pow((108*pow(c0,3)*pow(lz,6) + pow(tmax,3) 
      + sqrt(-pow(tmax,6) + pow((108*pow(c0,3)*pow(lz,6) + pow(tmax,3)),2))),1./3.)/(6*c0*MathUtility::sqr(lz));
    
    double s11FromEquilibrium = tmax/(lambdaValue*lz*lz);
    
    double optimalP = 2./3*c0*(pow(lambdaValue,3)-1)/lambdaValue;
    
    // constant load per surface area, t = that, but total load varying
    double that = 0.5;
    lambdaValue = 0.30285343213869*(1.81712059283214*that + pow((18.0*pow(c0,1.5) + 2.44948974278318*pow((54.0*pow(c0,3) - pow(that,3)),0.5)),2./3.))
      /(sqrt(c0)*pow((18.0*pow(c0,1.5) + 2.44948974278318*pow((54.0*pow(c0,3) - pow(that,3)),0.5)),1./3.));
    
    //double s11FromEquilibrium = 2*c0 + 4*c1/lambdaValue;
    
    //lambdaValue = 1.8;
    
    LOG(DEBUG) << "c0: " << c0 << ", c1: " << c1 << ", tmax: " << tmax ;
    LOG(DEBUG) << "lambdaValue: " << lambdaValue << ", lambdaValue^2: " << MathUtility::sqr(lambdaValue) << ", lambdaValue^{-1}: " << 1./lambdaValue << ", lambdaValue^{-2}: " << 1./MathUtility::sqr(lambdaValue);
    LOG(DEBUG) << "solution p value: " << optimalP;
    
    LOG(DEBUG) << "expected F: diag( " << lambdaValue << ", " << 1./sqrt(lambdaValue) << ", " << 1./sqrt(lambdaValue) << ")";
    LOG(DEBUG) << "expected C: diag( " << lambdaValue*lambdaValue << ", " << 1./lambdaValue << ", " << 1./lambdaValue << ")";
    LOG(DEBUG) << "expected I1: " << lambdaValue*lambdaValue + 2/lambdaValue;
    LOG(DEBUG) << "expected Sbar_00: " << 2*c0 + 4*c1/lambdaValue << ", C:Sbar = " << 2*c0*lambdaValue*lambdaValue + 4*c0/lambdaValue;
    LOG(DEBUG) << "expected correct S: diag( " << 2*c0 + 4*c1/lambdaValue << ", " << -2*c1/lambdaValue << ", " << -2*c1/lambdaValue <<  " ), "
     << " S penalty: diag( " << 4./3 * c0 * (1 - pow(lambdaValue,-3)) << ", " << 2./3*c0*(1 - pow(lambdaValue,3)) << ")";
    
    
    LOG(DEBUG) << "expected W[9]: " << s11FromEquilibrium * lambdaValue * MathUtility::sqr(lz) / 4.;
    
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
     
    if (this->data_.computeWithReducedVectors())
    {
      reduceVector(displacements, this->data_.displacementsReduced());
    }
  }
  else if(BasisOnMeshType::dim() == 2)
  {
    
    double lx = 1.0;
    double ly = 1.0;
    double tmax = 0.5;
    
    const double c0 = SEMT::Parameter<0>::get_value();
    const double c1 = SEMT::Parameter<1>::get_value();
    assert(c1 == 0);
    
    double lambdaValue = 1.0;
    
    // constant total load
    
    // analytical solution from equilibrium
    lambdaValue = (0.125*tmax + 0.0416666666666667*sqrt(34.6139896873778*pow(c0,1.33333333333333)*pow(ly,1.33333333333333)
    *pow((-9.0*pow(tmax,2.) + 1.73205080756888*pow((4096.0*pow(c0,4.)*pow(ly,4.) + 27.0*pow(tmax,4.)),0.5)),0.333333333333333)
    - 798.752188051931*pow(c0,2.66666666666667)*pow(ly,2.66666666666667)*pow((-9.0*pow(tmax,2.) 
    + 1.73205080756888*pow((4096.0*pow(c0,4.)*pow(ly,4.) + 27.0*pow(tmax,4.)),0.5)),-0.333333333333333) 
    + 9.0*pow(tmax,2.)) + 0.0416666666666667*sqrt(-34.6139896873778*pow(c0,1.33333333333333)*pow(ly,1.33333333333333)
    *pow((-9.0*pow(tmax,2.) + 1.73205080756888*pow((4096.0*pow(c0,4.)*pow(ly,4.) + 27.0*pow(tmax,4.)),0.5)),0.333333333333333)
    + 798.752188051931*pow(c0,2.66666666666667)*pow(ly,2.66666666666667)*pow((-9.0*pow(tmax,2.) 
    + 1.73205080756888*pow((4096.0*pow(c0,4.)*pow(ly,4.) + 27.0*pow(tmax,4.)),0.5)),-0.333333333333333)
    + 31.1769145362398*pow(tmax,3.)*pow((11.5379965624593*pow(c0,1.33333333333333)*pow(ly,1.33333333333333)*pow((-9.0*pow(tmax,2.) 
    + 1.73205080756888*pow((4096.0*pow(c0,4.)*pow(ly,4.) + 27.0*pow(tmax,4.)),0.5)),0.333333333333333) 
    - 266.250729350644*pow(c0,2.66666666666667)*pow(ly,2.66666666666667)*pow((-9.0*pow(tmax,2.)
    + 1.73205080756888*pow((4096.0*pow(c0,4.)*pow(ly,4.) + 27.0*pow(tmax,4.)),0.5)),-0.333333333333333)
    + 3.0*pow(tmax,2.)),-0.5) + 18.0*pow(tmax,2.)))/(c0*ly);
    
    // solution from penalty formulation (same)
    
    double optimalP = c0*(pow(lambdaValue,2) - pow(lambdaValue,-2.));
    
    LOG(DEBUG) << "c0: " << c0 << ", c1: " << c1 << ", tmax: " << tmax ;
    LOG(DEBUG) << "lambdaValue: " << lambdaValue << ", lambdaValue^2: " << MathUtility::sqr(lambdaValue) << ", lambdaValue^{-1}: " << 1./lambdaValue << ", lambdaValue^{-2}: " << 1./MathUtility::sqr(lambdaValue);
    LOG(DEBUG) << "solution p value: " << optimalP;
    
    std::vector<double> knownValues = {
     0.0, 0.0, 0.0, 
     (lambdaValue-1.0)*lx, 0.0, 0.0,
     0.0, (1./lambdaValue - 1.0)*ly, 0.0, 
     (lambdaValue-1.0)*lx, (1./lambdaValue - 1.0)*ly, 0.0
    };
    PetscUtility::setVector(knownValues, displacements);
    
    if (this->data_.computeWithReducedVectors())
    {
      reduceVector(displacements, this->data_.displacementsReduced());
    }
  }
  // set prescribed Dirchlet BC displacements values
  //applyDirichletBoundaryConditionsInDisplacements();
  
#if 0  
  VLOG(1) << "displacements values: " << PetscUtility::getStringVector(displacements);
  
  this->computeExternalVirtualWork(externalVirtualWork);
  VLOG(1) << "W_ext: " << PetscUtility::getStringVector(externalVirtualWork);
  
  
  
  this->computeInternalVirtualWork(internalVirtualWork);
  VLOG(1) << "W_int: " << PetscUtility::getStringVector(internalVirtualWork);
  
  std::vector<double> wExt, wInt;
  PetscUtility::getVectorEntries(externalVirtualWork, wExt);
  PetscUtility::getVectorEntries(internalVirtualWork, wInt);
  
  std::stringstream s;
  s << std::endl << "no.  Wext  Wint" << std::endl;
  for (int i=0; i<wExt.size(); i++)
  {
    s << i << "   " << wExt[i] << "   " << wInt[i] << std::endl;
  }
  LOG(DEBUG) << s.str();
#endif  
  LOG(DEBUG) << "------- debug method end, exit------";
  //exit(0);
  
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
  
  //debug();
  
  bool useAnalyticJacobian = PythonUtility::getOptionBool(this->specificSettings_, "analyticJacobian", true);
  PetscErrorCode ierr;
  
  Mat tangentStiffnessMatrix = this->data_.tangentStiffnessMatrix();
  Vec residual = this->data_.residual().values();
  Vec solverDisplacementVariable = this->data_.displacements().values();
  
  if (this->data_.computeWithReducedVectors())
  {
    residual = this->data_.residualReduced();
    solverDisplacementVariable = this->data_.displacementsReduced();
    tangentStiffnessMatrix = this->data_.tangentStiffnessMatrixReduced();
  }
  
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
  ierr = VecSet(solverDisplacementVariable, 0.0); CHKERRV(ierr);
  
  debug();
  
  // set prescribed Dirchlet BC displacements values
  applyDirichletBoundaryConditionsInDisplacements();
  
  VLOG(1) << "Dirichlet BC: indices: " << this->dirichletIndices_ << ", values: " << dirichletValues_;
  VLOG(1) << "initial values: " << PetscUtility::getStringVector(solverDisplacementVariable);

  if (!useAnalyticJacobian)
  {
     LOG(DEBUG) << "compute analytical jacobian to initialize non-zero structure of matrix";
     callbackJacobianAnalytic(*snes, solverDisplacementVariable, tangentStiffnessMatrix, tangentStiffnessMatrix, this);
  }
  
  LOG(DEBUG) << "------------------  start solve  ------------------";
  
  // solve the system nonlinearFunction(displacements) = 0
  // not sure if displacements has to be a different vector from the one used in the provided functions
  ierr = SNESSolve(*snes, NULL, solverDisplacementVariable); CHKERRV(ierr);
  
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
  
  // potentially expand solution vector
  if (this->data_.computeWithReducedVectors())
  {
    expandVector(solverDisplacementVariable, this->data_.displacements().values());
  }
  else 
  {
    assert(solverDisplacementVariable == this->data_.displacements().values());
  }
    
  // update geometry field from displacements, w = alpha*x+y
  VecWAXPY(this->data_.geometryActual().values(), 1.0, this->data_.geometryReference().values(), this->data_.displacements().values());
}
  
};