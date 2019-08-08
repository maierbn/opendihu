#include "specialized_solver/hyperelasticity/quasi_static_hyperelasticity_solver.h"

#include <Python.h>  // has to be the first included header

#include "specialized_solver/hyperelasticity/petsc_callbacks.h"
#include "solver/nonlinear.h"
#include "control/dihu_context.h"
#include "solver/solver_manager.h"

namespace TimeSteppingScheme
{

void QuasiStaticHyperelasticitySolver::
nonlinearSolve()
{
  LOG(TRACE) << "nonlinear solve";

  // create nonlinear solver PETSc context (snes)
  std::shared_ptr<Solver::Nonlinear> nonlinearSolver = this->context_.solverManager()->template solver<Solver::Nonlinear>(
    this->specificSettings_, this->displacementsFunctionSpace_->meshPartition()->mpiCommunicator());
  std::shared_ptr<SNES> snes = nonlinearSolver->snes();
  std::shared_ptr<KSP> ksp = nonlinearSolver->ksp();

  assert(snes != nullptr);

  this->initializeSolutionVariable();

  //VLOG(1) << "initial values: " << PetscUtility::getStringVector(solverVariableSolution);

  /*
  if (!useAnalyticJacobian_)
  {
     LOG(DEBUG) << "compute analytic jacobian to initialize non-zero structure of matrix";

     // set the displacement variables according to the values in the solve variable
     setFromSolverVariableSolution(solverVariableSolution);

     // compute tangent stiffness matrix for the first time. This constructs and initialized a Petsc Mat in solverMatrixTangentStiffness which will be used for all further computeJacobianAnalytic calls
     computeAnalyticStiffnessMatrix(solverMatrixTangentStiffness);

     LOG(DEBUG) << "done, now object is " << solverMatrixTangentStiffness;
  }
  */

  // set callback functions that are defined in petsc_callbacks.h
  typedef QuasiStaticHyperelasticitySolver ThisClass;

  PetscErrorCode (*callbackNonlinearFunction)(SNES, Vec, Vec, void *)              = *nonlinearFunction<ThisClass>;
  PetscErrorCode (*callbackJacobianAnalytic)(SNES, Vec, Mat, Mat, void *)          = *jacobianFunctionAnalytic<ThisClass>;
  PetscErrorCode (*callbackJacobianFiniteDifferences)(SNES, Vec, Mat, Mat, void *) = *jacobianFunctionFiniteDifferences<ThisClass>;
  PetscErrorCode (*callbackJacobianCombined)(SNES, Vec, Mat, Mat, void *)          = *jacobianFunctionCombined<ThisClass>;
  PetscErrorCode (*callbackMonitorFunction)(SNES, PetscInt, PetscReal, void *)     = *monitorFunction<ThisClass>;

  // set function
  PetscErrorCode ierr;
  ierr = SNESSetFunction(*snes, solverVariableResidual_, callbackNonlinearFunction, this); CHKERRV(ierr);

  // set jacobian
  if (useAnalyticJacobian_)
  {
    if (useNumericJacobian_)   // use combination of analytic jacobian also with finite differences
    {
      Mat solverMatrixTangentStiffnessFiniteDifferences;
      ierr = MatDuplicate(solverMatrixTangentStiffness_, MAT_DO_NOT_COPY_VALUES, &solverMatrixTangentStiffnessFiniteDifferences); CHKERRV(ierr);
      ierr = SNESSetJacobian(*snes, solverMatrixTangentStiffnessFiniteDifferences, solverMatrixTangentStiffness_, callbackJacobianCombined, this); CHKERRV(ierr);
      LOG(DEBUG) << "Use combination of numeric and analytic jacobian: " << solverMatrixTangentStiffness_;
    }
    else    // use pure analytic jacobian, without fd
    {
      ierr = SNESSetJacobian(*snes, solverMatrixTangentStiffness_, solverMatrixTangentStiffness_, callbackJacobianAnalytic, this); CHKERRV(ierr);
      LOG(DEBUG) << "Use only analytic jacobian: " << solverMatrixTangentStiffness_;
    }
  }
  else
  {
    // set function to compute jacobian from finite differences
    ierr = SNESSetJacobian(*snes, solverMatrixTangentStiffness_, solverMatrixTangentStiffness_, callbackJacobianFiniteDifferences, this); CHKERRV(ierr);
    LOG(DEBUG) << "Use Finite-Differences approximation for jacobian";
  }

  // prepare log file
  std::shared_ptr<std::ofstream> logFile = nullptr;
  if (this->specificSettings_.hasKey("residualNormLogFilename"))
  {
    std::string logFileName = this->specificSettings_.getOptionString("residualNormLogFilename", "residual_norm.txt");

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
  ierr = SNESSolve(*snes, NULL, solverVariableSolution_); CHKERRV(ierr);

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

  LOG(DEBUG) << "solution: " << getString(solverVariableSolution_);
  checkSolution(solverVariableSolution_);

  // set the displacement variables according to the values in the solve variable
  //setFromSolverVariableSolution(solverVariableSolution);
}

void QuasiStaticHyperelasticitySolver::
initializeSolutionVariable()
{
  if (this->useNestedMat_)
  {
    // zero initial values
    this->data_.displacements()->zeroEntries();
    this->data_.pressure()->zeroEntries();

    // set prescribed Dirchlet BC displacements values
    this->dirichletBoundaryConditions_->applyInVector(this->data_.displacements());
  }
  else
  {
    VecZeroEntries(solverVariableSolution_);

    // set x to x0 for values with dirichlet boundary condition
    PetscErrorCode ierr;

    typedef SpatialDiscretization::BoundaryConditionsBase<DisplacementsFunctionSpace,3>::BoundaryConditionsForComponent BoundaryConditionsForComponent;
    const std::array<BoundaryConditionsForComponent, 3> &boundaryConditionsByComponent = dirichletBoundaryConditions_->boundaryConditionsByComponent();

    for (int componentNo = 0; componentNo < 3; componentNo++)
    {
      const int nValues = boundaryConditionsByComponent[componentNo].dofNosLocal.size();

      std::vector<int> dofNos(nValues);
      std::vector<double> values(nValues);

      for (int i = 0; i < nValues; i++)
      {
        dofNos[i] = componentNo * displacementsFunctionSpace_->nDofsGlobal() + boundaryConditionsByComponent[componentNo].dofNosLocal[i];
      }

      for (int i = 0; i < nValues; i++)
      {
        values[i] = boundaryConditionsByComponent[componentNo].values[i];
      }

      ierr = VecSetValues(solverVariableSolution_, nValues, dofNos.data(), values.data(), INSERT_VALUES); CHKERRV(ierr);
    }
  }

  VecAssemblyBegin(solverVariableSolution_);
  VecAssemblyEnd(solverVariableSolution_);

  LOG(DEBUG) << "after initialization: " << this->getString(solverVariableSolution_);
}

void QuasiStaticHyperelasticitySolver::
applyDirichletBoundaryConditionsInJacobian(Vec x, Mat jac)
{
  // zero rows and columns in jac for which dirichlet values are set, set diagonal to 1
  PetscErrorCode ierr;

  typedef SpatialDiscretization::BoundaryConditionsBase<DisplacementsFunctionSpace,3>::BoundaryConditionsForComponent BoundaryConditionsForComponent;

  const std::array<BoundaryConditionsForComponent, 3> &boundaryConditionsByComponent = dirichletBoundaryConditions_->boundaryConditionsByComponent();

  for (int componentNo = 0; componentNo < 3; componentNo++)
  {
    const int nValues = boundaryConditionsByComponent[componentNo].dofNosLocal.size();

    std::vector<int> dofNos(nValues);

    for (int i = 0; i < nValues; i++)
    {
      dofNos[i] = componentNo * displacementsFunctionSpace_->nDofsGlobal() + boundaryConditionsByComponent[componentNo].dofNosLocal[i];
    }

    ierr = MatZeroRowsColumns(jac, nValues, dofNos.data(), 1.0, NULL, NULL); CHKERRV(ierr);
  }

  VLOG(2) << "jac with BC: " << PetscUtility::getStringMatrix(jac);
}

void QuasiStaticHyperelasticitySolver::
applyDirichletBoundaryConditionsInNonlinearFunction(Vec x, Vec f)
{
  // set f to x-x0 for values with dirichlet boundary condition
  PetscErrorCode ierr;

  typedef SpatialDiscretization::BoundaryConditionsBase<DisplacementsFunctionSpace,3>::BoundaryConditionsForComponent BoundaryConditionsForComponent;

  const std::array<BoundaryConditionsForComponent, 3> &boundaryConditionsByComponent = dirichletBoundaryConditions_->boundaryConditionsByComponent();

  for (int componentNo = 0; componentNo < 3; componentNo++)
  {
    const int nValues = boundaryConditionsByComponent[componentNo].dofNosLocal.size();

    std::vector<int> dofNos(nValues);
    std::vector<double> values(nValues), xValues(nValues);

    for (int i = 0; i < nValues; i++)
    {
      dofNos[i] = componentNo * displacementsFunctionSpace_->nDofsGlobal() + boundaryConditionsByComponent[componentNo].dofNosLocal[i];
    }

    ierr = VecGetValues(x, nValues, dofNos.data(), xValues.data()); CHKERRV(ierr);

    for (int i = 0; i < nValues; i++)
    {
      values[i] = xValues[i] - boundaryConditionsByComponent[componentNo].values[i];
      VLOG(3) << "  " << std::string(1,('x' + (char)componentNo)) << i << " set " << xValues[i] << "-" << boundaryConditionsByComponent[componentNo].values[i] << " = " << values[i];
    }

    ierr = VecSetValues(f, nValues, dofNos.data(), values.data(), INSERT_VALUES); CHKERRV(ierr);
  }
  VecAssemblyBegin(f);
  VecAssemblyEnd(f);

  VLOG(1) << "f with BC: " << this->getString(f);
}

void QuasiStaticHyperelasticitySolver::
evaluateNonlinearFunction(Vec x, Vec f)
{
  VLOG(1) << "evaluateNonlinearFunction at " << getString(x);

  // VecAXPY(y, alpha, x), y = alpha x + y.
  VecCopy(x, f);
  VecPointwiseMult(f,f,f);   // VecPointwiseMult(w,x,y), w = x*y
  VecScale(f, 2.0);
  VecShift(f, -0.5);
  // f(x) = 2*x^2 - 0.5   =>   f'(x) = 2*x,  f(0.5) = 0
}

void QuasiStaticHyperelasticitySolver::
evaluateAnalyticJacobian(Vec x, Mat jac)
{
  int nValues = 3*displacementsFunctionSpace_->nDofsGlobal() + pressureFunctionSpace_->nDofsGlobal();

  // set diagonal entries
  for (int i = 0; i < nValues; i++)
  {
    MatSetValue(jac, i, i, 3.0, INSERT_VALUES);
  }

  MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY);
}

void QuasiStaticHyperelasticitySolver::
checkSolution(Vec x)
{
  if (this->useNestedMat_)
  {
  }
  else
  {
    // check if Dirichlet BC are met

    typedef SpatialDiscretization::BoundaryConditionsBase<DisplacementsFunctionSpace,3>::BoundaryConditionsForComponent BoundaryConditionsForComponent;
    const std::array<BoundaryConditionsForComponent, 3> &boundaryConditionsByComponent = dirichletBoundaryConditions_->boundaryConditionsByComponent();

    double l2NormDirichletBC = 0.0;
    int nEntries = 0;
    for (int componentNo = 0; componentNo < 3; componentNo++)
    {
      const int nValues = boundaryConditionsByComponent[componentNo].dofNosLocal.size();

      std::vector<int> dofNos(nValues);
      std::vector<double> values(nValues), xValues(nValues);

      for (int i = 0; i < nValues; i++)
      {
        dofNos[i] = componentNo * displacementsFunctionSpace_->nDofsGlobal() + boundaryConditionsByComponent[componentNo].dofNosLocal[i];
      }

      PetscErrorCode ierr;
      ierr = VecGetValues(x, nValues, dofNos.data(), xValues.data()); CHKERRV(ierr);

      for (int i = 0; i < nValues; i++)
      {
        double contribution = MathUtility::sqr(xValues[i] - boundaryConditionsByComponent[componentNo].values[i]);
        l2NormDirichletBC += contribution;
        nEntries++;

        if (fabs(contribution) > 1e-5)
          LOG(DEBUG) << "Component " << componentNo << " dof " << i << ", value is " << xValues[i] << ", should be " << boundaryConditionsByComponent[componentNo].values[i];
      }
    }

    l2NormDirichletBC /= nEntries;
    LOG(DEBUG) << "L2-norm Dirchlet-BC residual: " << l2NormDirichletBC;

    if (l2NormDirichletBC < 1e-5)
    {
      LOG(DEBUG) << "Dirichlet BC are met.";
    }
    else
    {
      LOG(DEBUG) << "Dirichlet BC are NOT met.";
    }

    // check if function is zero
    evaluateNonlinearFunction(x, solverVariableResidual_);

    // set rows in f to x - x0 for which dirichlet BC in x is given
    applyDirichletBoundaryConditionsInNonlinearFunction(x, solverVariableResidual_);

    int nValues = 3*displacementsFunctionSpace_->nDofsGlobal() + pressureFunctionSpace_->nDofsGlobal();
    for (int i = 0; i < nValues; i++)
    {
      double value;
      PetscErrorCode ierr;
      ierr = VecGetValues(solverVariableResidual_, 1, &i, &value); CHKERRV(ierr);
      if (fabs(value) > 1e-5)
        LOG(DEBUG) << "Dof " << i << ", residual is " << value << " should be zero.";
    }

    PetscReal l2NormResidual;
    VecNorm(solverVariableResidual_, NORM_2, &l2NormResidual);
    LOG(DEBUG) << "L2-norm residual: " << l2NormResidual;

    if (l2NormResidual < 1e-5)
    {
      LOG(DEBUG) << "Root found.";
    }
    else
    {
      LOG(DEBUG) << "Root NOT found.";
    }
  }
}


std::string QuasiStaticHyperelasticitySolver::
getString(Vec x)
{
  if (this->useNestedMat_)
  {
  }
  else
  {
    std::stringstream result;

    std::array<std::vector<double>,4> values;
    for (int componentNo = 0; componentNo < 4; componentNo++)
    {
      int nValues = displacementsFunctionSpace_->nDofsGlobal();
      if (componentNo == 3)
        nValues = pressureFunctionSpace_->nDofsGlobal();

      values[componentNo].resize(nValues);

      std::vector<int> dofNos(nValues);

      for (int i = 0; i < nValues; i++)
      {
        dofNos[i] = componentNo * displacementsFunctionSpace_->nDofsGlobal() + i;
      }

      PetscErrorCode ierr;
      ierr = VecGetValues(x, nValues, dofNos.data(), values[componentNo].data()); CHKERRABORT(displacementsFunctionSpace_->meshPartition()->mpiCommunicator(), ierr);
    }

    result << "[x:" << values[0] << " y:" << values[1] << " z:" << values[2] << " p:" << values[3] << "]";
    return result.str();
  }
  return std::string("no");
}

} // namespace TimeSteppingScheme
