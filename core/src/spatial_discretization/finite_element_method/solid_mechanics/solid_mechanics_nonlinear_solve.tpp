#include "spatial_discretization/finite_element_method/solid_mechanics/02_stiffness_matrix_incompressible.h"

#include <Python.h>  // this has to be the first included header
#include <iostream>
#include <petscsnes.h>
#include <petscksp.h>

#include "easylogging++.h"
#include "solver/nonlinear.h"
#include "data_management/finite_elements_solid_mechanics.h"
#include "control/performance_measurement.h"

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

  // copy entries of the solution vector to displacements or displacements and pressure
  object->setFromSolverVariableSolution(u);

  object->applyDirichletBoundaryConditionsInDisplacements(object->data());

  // compute the lhs which is the virtual work
  object->evaluateNonlinearFunction(f);

  // set rows in f to 0 for which dirichlet BC in u is given
  object->applyDirichletBoundaryConditionsInNonlinearFunction(f, object->data());

  LOG(DEBUG) << "--     displ/pressure input with BC: " << PetscUtility::getStringVector(u);
  LOG(DEBUG) << "--         nonlinear output with BC: " << PetscUtility::getStringVector(f);

  // compute and output function norm
  PetscReal functionNorm;
  VecNorm(f, NORM_2, &functionNorm);
  LOG(DEBUG) << "function norm: " << functionNorm;

  //object->writeOutput();

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

  object->setFromSolverVariableSolution(x);

  LOG(DEBUG) << "tangentStiffnessMatrix jac: " << jac << ", b: " << b;

  // compute the tangent stiffness matrix
  object->computeAnalyticStiffnessMatrix(jac);

  //MatCopy(tangentStiffnessMatrix,jac,SAME_NONZERO_PATTERN);
  //MatCopy(tangentStiffnessMatrix,b,SAME_NONZERO_PATTERN);

  VLOG(1) << "-- computed tangent stiffness matrix analytically: " << PetscUtility::getStringMatrix(jac);
  VLOG(1) << "-- non-zeros pattern: " << std::endl << PetscUtility::getStringSparsityPattern(jac);

  object->writeOutput();
  return 0;
};


template<typename T>
PetscErrorCode jacobianFunctionFiniteDifferences(SNES snes, Vec x, Mat jac, Mat b, void *context)
{
  T* object = static_cast<T*>(context);

  SNESComputeJacobianDefault(snes, x, jac, b, context);

  // zero rows and columns for which Dirichlet BC is set
  object->applyDirichletBoundaryConditionsInStiffnessMatrix(b, object->data());
  VLOG(1) << "-- computed tangent stiffness matrix by finite differences: " << PetscUtility::getStringMatrix(jac);
  VLOG(1) << "-- non-zeros pattern: " << std::endl << PetscUtility::getStringSparsityPattern(jac);

  object->writeOutput();
  return 0;
};


template<typename T>
PetscErrorCode jacobianFunctionCombined(SNES snes, Vec x, Mat jac, Mat b, void *context)
{
  T* object = static_cast<T*>(context);

  object->setFromSolverVariableSolution(x);

  // compute the finite differences jacobian in the main jacobian slot jac
  SNESComputeJacobianDefault(snes, x, jac, jac, context);

  // zero rows and columns for which Dirichlet BC is set
  object->applyDirichletBoundaryConditionsInStiffnessMatrix(jac, object->data());

  // compute the tangent stiffness matrix, stored in the preconditioner slot b
  object->computeAnalyticStiffnessMatrix(b);

  object->writeOutput();

  LOG(INFO) << "terminate in jacobianFunctionCombined";
  exit(0);
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

  // if log file was given, write residual norm to log file
  if (mctx != NULL)
  {
    std::ofstream *logFile = (std::ofstream *)(mctx);
    *logFile << its << ";" << norm << std::endl;
  }

  return 0;
};

template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void SolidMechanicsNonlinearSolve<BasisOnMeshType,QuadratureType,Term>::
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
      const int D = BasisOnMeshType::dim();
      const int nUnknowns = this->data_.mesh()->nDofs() * D;

      this->reduceVector(displacements, this->data_.solverVariableSolution(), nUnknowns);
    }

  }
  else if(BasisOnMeshType::dim() == 2)
  {

#if 0   // 1 element
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
     0.0, 0.0,
     (lambdaValue-1.0)*lx, 0.0,
     0.0, (1./lambdaValue - 1.0)*ly,
     (lambdaValue-1.0)*lx, (1./lambdaValue - 1.0)*ly
    };
    PetscUtility::setVector(knownValues, displacements);

    if (this->data_.computeWithReducedVectors())
    {
      const int D = BasisOnMeshType::dim();
      const int nUnknowns = this->data_.mesh()->nDofs() * D;

      this->reduceVector(displacements, this->data_.solverVariableSolution(), nUnknowns);
    }
  }
  // set prescribed Dirchlet BC displacements values
  //applyDirichletBoundaryConditionsInDisplacements(this->data_.computeWithReducedVectors());
#endif

#if 1   // 1x2 elements
   double analytic_lambda = 1.96511022361548;
   double lx_settings = 1.5;
   double ly_settings = 0.6;
   const double c0 = SEMT::Parameter<0>::get_value();

   double analyticPressure = c0 * (MathUtility::sqr(analytic_lambda) - 1./MathUtility::sqr(analytic_lambda));
   //double analyticPressure = 2*c0 / MathUtility::sqr(analytic_lambda);

   VLOG(1) << "analyticPressure: " << analyticPressure;

   std::vector<double> knownValues = {
     0.00,                             0.0,
     0.25*analytic_lambda*lx_settings, 0.0,
     0.50*analytic_lambda*lx_settings, 0.0,
     0.75*analytic_lambda*lx_settings, 0.0,
     1.00*analytic_lambda*lx_settings, 0.0,

     0.00,                             0.5/analytic_lambda*ly_settings,
     0.25*analytic_lambda*lx_settings, 0.5/analytic_lambda*ly_settings,
     0.50*analytic_lambda*lx_settings, 0.5/analytic_lambda*ly_settings,
     0.70*analytic_lambda*lx_settings, 0.5/analytic_lambda*ly_settings,
     1.00*analytic_lambda*lx_settings, 0.5/analytic_lambda*ly_settings,

     0.00,                             1.00/analytic_lambda*ly_settings,
     0.25*analytic_lambda*lx_settings, 1.00/analytic_lambda*ly_settings,
     0.50*analytic_lambda*lx_settings, 1.00/analytic_lambda*ly_settings,
     0.70*analytic_lambda*lx_settings, 1.00/analytic_lambda*ly_settings,
     1.00*analytic_lambda*lx_settings, 1.00/analytic_lambda*ly_settings,

     analyticPressure, analyticPressure, analyticPressure,
     analyticPressure, analyticPressure, analyticPressure
   };

   if (this->data_.computeWithReducedVectors())
   {
     const int nUnknowns = this->nUnknowns();

     std::vector<double> reducedVector(nUnknowns - this->dirichletIndices_.size());
     dof_no_t reducedIndex = 0;
     std::vector<dof_no_t>::const_iterator dirichletIndicesIter = this->dirichletIndices_.begin();

     for (dof_no_t currentDofNo = 0; currentDofNo < nUnknowns; currentDofNo++)
     {
       // exclude variables for which Dirichlet BC are set
       if (dirichletIndicesIter != this->dirichletIndices_.end())
       {
         if (currentDofNo == *dirichletIndicesIter)
         {
           dirichletIndicesIter++;
           continue;
         }
       }

       reducedVector[reducedIndex++] = knownValues[currentDofNo];
     }

     PetscUtility::setVector(reducedVector, this->data_.solverVariableSolution());
     this->setFromSolverVariableSolution(this->data_.solverVariableSolution());


     Vec solverVariableResidual = this->data_.solverVariableResidual();


     this->evaluateNonlinearFunction(solverVariableResidual);

     exit(0);
   }
#endif
  }
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

// general implementation of nonlinear solving for solid mechanics
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void SolidMechanicsNonlinearSolve<BasisOnMeshType,QuadratureType,Term>::
solve()
{
  LOG(TRACE) << "FiniteElementMethod::solve (nonlinear)";

  bool useAnalyticJacobian = PythonUtility::getOptionBool(this->specificSettings_, "analyticJacobian", true);
  bool useNumericJacobian = PythonUtility::getOptionBool(this->specificSettings_, "numericJacobian", true);

  if (!useAnalyticJacobian && !useNumericJacobian)
  {
    LOG(WARNING) << "Can not set both \"analyticJacobian\" and \"numericJacobian\" to False, using numeric jacobian.";
    useNumericJacobian = true;
  }

  PetscErrorCode ierr;

  Mat solverMatrixTangentStiffness = this->data_.solverMatrixTangentStiffness();
  Vec solverVariableResidual = this->data_.solverVariableResidual();
  Vec solverVariableSolution = this->data_.solverVariableSolution();

  if (!useAnalyticJacobian && !this->data_.computeWithReducedVectors())
  {
    LOG(FATAL) << "This combination of numeric jacobian and non-reduced vectors does not work.";
  }

  // create nonlinear solver PETSc context (snes)
  std::shared_ptr<Solver::Nonlinear> nonlinearSolver = this->context_.solverManager()->template solver<Solver::Nonlinear>(this->specificSettings_);
  std::shared_ptr<SNES> snes = nonlinearSolver->snes();
  std::shared_ptr<KSP> ksp = nonlinearSolver->ksp();

  assert(snes != nullptr);

  // zero initial values
  ierr = VecSet(solverVariableSolution, 0.0); CHKERRV(ierr);

  //debug();

  // set prescribed Dirchlet BC displacements values
  this->applyDirichletBoundaryConditionsInDisplacements(this->data_);

  VLOG(1) << "Dirichlet BC: indices: " << this->dirichletIndices_ << ", values: " << this->dirichletValues_;
  VLOG(1) << "initial values: " << PetscUtility::getStringVector(solverVariableSolution);

  if (!useAnalyticJacobian || (useAnalyticJacobian && this->data_.computeWithReducedVectors()))
  {
     LOG(DEBUG) << "compute analytic jacobian to initialize non-zero structure of matrix";

     // set the displacement variables according to the values in the solve variable
     setFromSolverVariableSolution(solverVariableSolution);

     // compute tangent stiffness matrix for the first time. This constructs and initialized a Petsc Mat in solverMatrixTangentStiffness which will be used for all further computeJacobianAnalytic calls
     computeAnalyticStiffnessMatrix(solverMatrixTangentStiffness);

     LOG(DEBUG) << "done, now object is " << solverMatrixTangentStiffness;
  }

  typedef FiniteElementMethodStiffnessMatrix<
    BasisOnMeshType,
    QuadratureType,
    Term,
    typename BasisOnMeshType::Mesh,
    Term,
    typename BasisOnMeshType::BasisFunction
  > ThisClass;

  PetscErrorCode (*callbackNonlinearFunction)(SNES, Vec, Vec, void *) = *nonlinearFunction<ThisClass>;
  PetscErrorCode (*callbackJacobianAnalytic)(SNES, Vec, Mat, Mat, void *) = *jacobianFunctionAnalytic<ThisClass>;
  PetscErrorCode (*callbackJacobianFiniteDifferences)(SNES, Vec, Mat, Mat, void *) = *jacobianFunctionFiniteDifferences<ThisClass>;
  PetscErrorCode (*callbackJacobianCombined)(SNES, Vec, Mat, Mat, void *) = *jacobianFunctionCombined<ThisClass>;
  PetscErrorCode (*callbackMonitorFunction)(SNES, PetscInt, PetscReal, void *) = *monitorFunction<ThisClass>;

  // set function
  ierr = SNESSetFunction(*snes, solverVariableResidual, callbackNonlinearFunction, this); CHKERRV(ierr);

  // set jacobian
  if (useAnalyticJacobian)
  {
    if(useNumericJacobian)   // use combination of analytic jacobian also with finite differences
    {
      Mat &solverMatrixTangentStiffnessFiniteDifferences = this->data_.solverMatrixTangentStiffnessFiniteDifferences();
      ierr = MatDuplicate(solverMatrixTangentStiffness, MAT_DO_NOT_COPY_VALUES, &solverMatrixTangentStiffnessFiniteDifferences); CHKERRV(ierr);
      ierr = SNESSetJacobian(*snes, solverMatrixTangentStiffnessFiniteDifferences, solverMatrixTangentStiffness, callbackJacobianCombined, this); CHKERRV(ierr);
      LOG(DEBUG) << "Use combination of numeric and analytic jacobian: " << solverMatrixTangentStiffness;
    }
    else    // use pure analytic jacobian, without fd
    {
      ierr = SNESSetJacobian(*snes, solverMatrixTangentStiffness, solverMatrixTangentStiffness, callbackJacobianAnalytic, this); CHKERRV(ierr);
      LOG(DEBUG) << "Use only analytic jacobian: " << solverMatrixTangentStiffness;
    }
  }
  else
  {
    // set function to compute jacobian from finite differences
    ierr = SNESSetJacobian(*snes, solverMatrixTangentStiffness, solverMatrixTangentStiffness, callbackJacobianFiniteDifferences, this); CHKERRV(ierr);
    LOG(DEBUG) << "Use Finite-Differences approximation for jacobian";
  }

  // prepare log file
  std::shared_ptr<std::ofstream> logFile = nullptr;
  if (PythonUtility::hasKey(this->specificSettings_, "logfile"))
  {
    std::string logFileName = PythonUtility::getOptionString(this->specificSettings_, "logfile", "residual_norm.txt");

    logFile = std::make_shared<std::ofstream>(logFileName, std::ios::out | std::ios::binary | std::ios::trunc);

    if (!logFile->is_open())
    {
      LOG(WARNING) << "Could not open log file for residual norm, \"" << logFileName << "\".";
      logFile = nullptr;
    }
  }

  // set monitor function
  ierr = SNESMonitorSet(*snes, callbackMonitorFunction, logFile.get(), NULL); CHKERRV(ierr);

  LOG(DEBUG) << "------------------  start solve  ------------------";

  // solve the system nonlinearFunction(displacements) = 0
  // not sure if displacements has to be a different vector from the one used in the provided functions
  ierr = SNESSolve(*snes, NULL, solverVariableSolution); CHKERRV(ierr);

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

  // close log file
  if (logFile != nullptr)
  {
    logFile->close();
  }

  // set the displacement variables according to the values in the solve variable
  setFromSolverVariableSolution(solverVariableSolution);

  updateGeometryActual();
}

// general implementation of nonlinear solving for solid mechanics
template<typename BasisOnMeshType, typename QuadratureType, typename Term>
void SolidMechanicsNonlinearSolve<BasisOnMeshType,QuadratureType,Term>::
updateGeometryActual()
{
  // update geometry field from displacements
  if (BasisOnMeshType::dim() == 2)  // 2D problem
  {
    const int nUnknowns3D = this->data_.mesh()->nDofs() * 3;

    // expand 2D vector to 3D vector fullIncrement
    this->expandVectorTo3D(this->data_.displacements().values(), this->data_.fullIncrement(), nUnknowns3D);

    // w = alpha*x+y
    VecWAXPY(this->data_.geometryActual().values(), 1.0, this->data_.geometryReference().values(), this->data_.fullIncrement());
  }
  else  // 3D problem
  {
    // w = alpha*x+y
    VecWAXPY(this->data_.geometryActual().values(), 1.0, this->data_.geometryReference().values(), this->data_.displacements().values());
  }
}

};