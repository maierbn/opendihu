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

  //VLOG(1) << "initial values: " << PetscUtility::getStringVectorVector(solverVariableSolution);

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

  nonlinearSolver->dumpMatrixRightHandSide(solverVariableResidual_);

  LOG(DEBUG) << "solution: " << getString(solverVariableSolution_);
  checkSolution(solverVariableSolution_);

  // set the displacement variables according to the values in the solve variable
  //setFromSolverVariableSolution(solverVariableSolution);
}

void QuasiStaticHyperelasticitySolver::
initializeSolutionVariable()
{
  // set variable to all zero and dirichlet boundary condition values
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
  if (this->useNestedMat_)
  {
    // already taken care of in construction of matrix

    // zeroRows would be possible,
    // zeroColumns is not possible to implement here, we need information
    // about the mesh.

    // this should not be needed, because the numeric jacobian has the dirichlet BC built-in because the nonlinear function has them

    // the analytic jacobian should consider Dirichlet BC within construction

    return;

    VLOG(2) << "before setting dirichlet BC jacobian:" << PetscUtility::getStringMatrix(jac);

    // for serial testing
    typedef SpatialDiscretization::BoundaryConditionsBase<DisplacementsFunctionSpace,3>::BoundaryConditionsForComponent BoundaryConditionsForComponent;

    const std::array<BoundaryConditionsForComponent, 3> &boundaryConditionsByComponent = dirichletBoundaryConditions_->boundaryConditionsByComponent();
    PetscErrorCode ierr;

    // allow new allocation of diagonal entries
    for (int i = 0; i < 4; i++)
    {
      for (int j = 0; j < 4; j++)
      {
        MatSetOption(submatrices_[i*4 + j],MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
      }
    }

    for (int componentNo = 0; componentNo < 3; componentNo++)
    {
      const int nValues = boundaryConditionsByComponent[componentNo].dofNosLocal.size();

      int ownershipBegin = 0;
      int ownershipEnd = 0;
      MatGetOwnershipRange(submatrices_[componentNo*4 + componentNo], &ownershipBegin, &ownershipEnd);

      std::vector<int> dofsNosGlobal(nValues);
      for (int i = 0; i < nValues; i++)
      {
        dofsNosGlobal[i] = ownershipBegin + boundaryConditionsByComponent[componentNo].dofNosLocal[i];
      }

      // own matrix
      ierr = MatZeroRowsColumns(submatrices_[componentNo*4 + componentNo], nValues, dofsNosGlobal.data(), 1.0, NULL, NULL); CHKERRV(ierr);

      MatAssemblyBegin(submatrices_[componentNo*4 + componentNo], MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(submatrices_[componentNo*4 + componentNo], MAT_FINAL_ASSEMBLY);

      // zero rows

      // loop over column blocks
      for (int j = 0; j < 4; j++)
      {
        if (j == componentNo)
          continue;

        ierr = MatZeroRows(submatrices_[componentNo*4 + j], nValues, dofsNosGlobal.data(), 0.0, NULL, NULL); CHKERRV(ierr);

        MatAssemblyBegin(submatrices_[componentNo*4 + j], MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(submatrices_[componentNo*4 + j], MAT_FINAL_ASSEMBLY);
      }

      // zero columns
      // loop over row blocks
      for (int i = 0; i < 4; i++)
      {
        if (i == componentNo)
          continue;

        int nRows = displacementsFunctionSpace_->nDofsGlobal();
        if (i == 3)
          nRows = pressureFunctionSpace_->nDofsGlobal();

        std::vector<double> values(nValues, 0.0);

        // loop over actual rows
        for (int j = 0; j < nRows; j++)
        {
          ierr = MatSetValues(submatrices_[i*4 + componentNo], 1, &j, nValues, dofsNosGlobal.data(), values.data(), INSERT_VALUES); CHKERRV(ierr);
        }

        MatAssemblyBegin(submatrices_[i*4 + componentNo], MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(submatrices_[i*4 + componentNo], MAT_FINAL_ASSEMBLY);
      }
    }

    VLOG(2) << "after setting dirichlet BC jacobian:" << PetscUtility::getStringMatrix(jac);
  }
  else
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
}

void QuasiStaticHyperelasticitySolver::
applyDirichletBoundaryConditionsInNonlinearFunction(Vec x, Vec f)
{
  // set f to x-x0 for values with dirichlet boundary condition

  typedef SpatialDiscretization::BoundaryConditionsBase<DisplacementsFunctionSpace,3>::BoundaryConditionsForComponent BoundaryConditionsForComponent;
  const std::array<BoundaryConditionsForComponent, 3> &boundaryConditionsByComponent = dirichletBoundaryConditions_->boundaryConditionsByComponent();

  PetscErrorCode ierr;

  if (this->useNestedMat_)
  {
    // get the subvectors of the nested Vec, the 4 subvectors correspond to [ux,uy,uz,p]
    Vec *xSubvectors;
    Vec *fSubvectors;

    int nSubvecs;
    VecNestGetSubVecs(x, &nSubvecs, &xSubvectors);
    VecNestGetSubVecs(f, &nSubvecs, &fSubvectors);

    // loop over the first 3 components (x,y,z) and set the value of the function to f(x) = x-x0
    for (int componentNo = 0; componentNo < 3; componentNo++)
    {
      const int nValues = boundaryConditionsByComponent[componentNo].dofNosLocal.size();
      std::vector<double> xValues(nValues);

      int ownershipBegin = 0;
      int ownershipEnd = 0;
      VecGetOwnershipRange(xSubvectors[componentNo], &ownershipBegin, &ownershipEnd);

      std::vector<int> dofsNosGlobal(nValues);
      for (int i = 0; i < nValues; i++)
      {
        dofsNosGlobal[i] = ownershipBegin + boundaryConditionsByComponent[componentNo].dofNosLocal[i];
      }

      // get the values x
      ierr = VecGetValues(xSubvectors[componentNo], nValues, dofsNosGlobal.data(), xValues.data());

      // compute the values f = x - x0
      std::vector<double> values(nValues);
      for (int i = 0; i < nValues; i++)
      {
        values[i] = xValues[i] - boundaryConditionsByComponent[componentNo].values[i];
      }

      // set the values f
      ierr = VecSetValues(fSubvectors[componentNo], nValues, dofsNosGlobal.data(), values.data(), INSERT_VALUES); CHKERRV(ierr);

      VecAssemblyBegin(fSubvectors[componentNo]);
      VecAssemblyEnd(fSubvectors[componentNo]);
    }
  }
  else
  {
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
  }

  VLOG(1) << "f with BC: " << this->getString(f);
}

void QuasiStaticHyperelasticitySolver::
evaluateNonlinearFunction(Vec x, Vec f)
{
  VLOG(1) << "evaluateNonlinearFunction at " << getString(x);

  if (this->useNestedMat_)
  {
    bool backupVecs = f != solverVariableResidual_ || x != solverVariableSolution_;

    // backup the values of solverVariableSolution_ and solverVariableResidual_ to be able to work with them now
    if (backupVecs)
    {
      // this happens in computation of the numeric jacobian
      VecSwap(x, solverVariableSolution_);
      VecSwap(f, solverVariableResidual_);
    }

    // evaluate nonlinear function F(X), where X=solverVariableSolution_ and F=solverVariableResidual_
    VecCopy(solverVariableSolution_, solverVariableResidual_);
    VecPointwiseMult(solverVariableResidual_,solverVariableResidual_,solverVariableResidual_);   // VecPointwiseMult(w,x,y), w = x*y
    VecScale(solverVariableResidual_, 2.0);
    VecShift(solverVariableResidual_, -0.5);

    // restore the values of solverVariableSolution_ and solverVariableResidual_ to their original values
    if (backupVecs)
    {
      VecSwap(x, solverVariableSolution_);
      VecSwap(f, solverVariableResidual_);
    }
  }
  else
  {
    // VecAXPY(y, alpha, x), y = alpha x + y.
    VecCopy(x, f);
    VecPointwiseMult(f,f,f);   // VecPointwiseMult(w,x,y), w = x*y
    VecScale(f, 2.0);
    VecShift(f, -0.5);
    // test 1: f(x) = 2*x - 1   =>   f'(x) = 2,  f(0.5) = 0
    // test 2: f(x) = 2*x^2 - 0.5   =>   f'(x) = 4*x,  f(0.5) = 0
  }
}

void QuasiStaticHyperelasticitySolver::
evaluateAnalyticJacobian(Vec x, Mat jac)
{
  VLOG(1) << "evaluateAnalyticJacobian";
  if (this->useNestedMat_)
  {
    typedef SpatialDiscretization::BoundaryConditionsBase<DisplacementsFunctionSpace,3>::BoundaryConditionsForComponent BoundaryConditionsForComponent;
    const std::array<BoundaryConditionsForComponent, 3> &boundaryConditionsByComponent = dirichletBoundaryConditions_->boundaryConditionsByComponent();

    for (int componentNo = 0; componentNo < 4; componentNo++)
    {
      int nValues = displacementsFunctionSpace_->nDofsGlobal();
      if (componentNo == 3)
        nValues = pressureFunctionSpace_->nDofsGlobal();

      // set diagonal entries
      int j = 0;
      for (int i = 0; i < nValues; i++)
      {
        // skip dofs with Dirichlet BC
        if (componentNo < 3)
        {
          if (j < boundaryConditionsByComponent[componentNo].dofNosLocal.size())
          {
            if (i == boundaryConditionsByComponent[componentNo].dofNosLocal[j])
            {
              j++;
              continue;
            }
          }
        }

        double value = 2.0;
        VLOG(1) << "matrix " << componentNo << "," << componentNo << " (" << componentNo*4 + componentNo << ") set entry " << i << "," << i << " to " << value;

        PetscErrorCode ierr;
        ierr = MatSetValue(submatrices_[componentNo*4 + componentNo], i, i, value, INSERT_VALUES); CHKERRV(ierr);
        //MatSetValuesLocal(submatrices_[componentNo*4 + componentNo], 1, &i, 1, &i, &value, INSERT_VALUES);
      }
    }
    MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY);

    VLOG(2) << "analytic jacobian with BC: " << PetscUtility::getStringMatrix(jac);
  }
  else
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
}

void QuasiStaticHyperelasticitySolver::
checkSolution(Vec x)
{
  if (this->useNestedMat_)
  {
    // check if Dirichlet BC are met

    typedef SpatialDiscretization::BoundaryConditionsBase<DisplacementsFunctionSpace,3>::BoundaryConditionsForComponent BoundaryConditionsForComponent;
    const std::array<BoundaryConditionsForComponent, 3> &boundaryConditionsByComponent = dirichletBoundaryConditions_->boundaryConditionsByComponent();

    double l2NormDirichletBC = 0.0;
    int nEntries = 0;
    for (int componentNo = 0; componentNo < 3; componentNo++)
    {
      const int nValues = boundaryConditionsByComponent[componentNo].dofNosLocal.size();

      std::vector<double> values(nValues), xValues(nValues);

      int ownershipBegin = 0;
      int ownershipEnd = 0;
      VecGetOwnershipRange(subvectorsSolution_[componentNo], &ownershipBegin, &ownershipEnd);

      std::vector<int> dofsNosGlobal(nValues);
      for (int i = 0; i < nValues; i++)
      {
        dofsNosGlobal[i] = ownershipBegin + boundaryConditionsByComponent[componentNo].dofNosLocal[i];
      }

      PetscErrorCode ierr;
      ierr = VecGetValues(subvectorsSolution_[componentNo], nValues, dofsNosGlobal.data(), xValues.data()); CHKERRV(ierr);

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
      LOG(DEBUG) << ANSI_COLOR_GREEN "Dirichlet BC are met." ANSI_COLOR_RESET;
    }
    else
    {
      LOG(DEBUG) << ANSI_COLOR_RED "Dirichlet BC are NOT met." ANSI_COLOR_RESET;
    }

    // check if function is zero
    evaluateNonlinearFunction(x, solverVariableResidual_);

    // set rows in f to x - x0 for which dirichlet BC in x is given
    applyDirichletBoundaryConditionsInNonlinearFunction(x, solverVariableResidual_);

    double l2NormResidualSum = 0;
    for (int componentNo = 0; componentNo < 4; componentNo++)
    {
      int nValues = displacementsFunctionSpace_->nDofsLocalWithoutGhosts();
      if (componentNo == 3)
        nValues = pressureFunctionSpace_->nDofsLocalWithoutGhosts();

      int ownershipBegin = 0;
      int ownershipEnd = 0;
      VecGetOwnershipRange(subvectorsResidual_[componentNo], &ownershipBegin, &ownershipEnd);

      for (int i = 0; i < nValues; i++)
      {
        if (ownershipBegin <= i && i < ownershipEnd)
        {
          double value;
          PetscErrorCode ierr;
          ierr = VecGetValues(subvectorsResidual_[componentNo], 1, &i, &value); CHKERRV(ierr);
          if (fabs(value) > 1e-5)
            LOG(DEBUG) << std::string(1,char('x'+componentNo)) << " " << i << ", residual is " << value << " should be zero.";
        }
      }

      PetscReal l2NormResidual;
      VecNorm(subvectorsResidual_[componentNo], NORM_2, &l2NormResidual);
      LOG(DEBUG) << "L2-norm residual: " << std::string(1,char('x'+componentNo)) << l2NormResidual;
      l2NormResidualSum += l2NormResidual;
    }

    if (l2NormResidualSum < 1e-5)
    {
      LOG(DEBUG) << ANSI_COLOR_GREEN "Root found." ANSI_COLOR_RESET;
    }
    else
    {
      LOG(DEBUG) << ANSI_COLOR_RED "Root NOT found." ANSI_COLOR_RESET;
    }
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
      LOG(DEBUG) << ANSI_COLOR_GREEN "Dirichlet BC are met." ANSI_COLOR_RESET;
    }
    else
    {
      LOG(DEBUG) << ANSI_COLOR_RED "Dirichlet BC are NOT met." ANSI_COLOR_RESET;
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
      LOG(DEBUG) << ANSI_COLOR_GREEN "Root found." ANSI_COLOR_RESET;
    }
    else
    {
      LOG(DEBUG) << ANSI_COLOR_RED "Root NOT found." ANSI_COLOR_RESET;
    }
  }
}


std::string QuasiStaticHyperelasticitySolver::
getString(Vec x)
{
  if (this->useNestedMat_)
  {
    std::stringstream result;
    result << "[xyz: " << *this->data_.displacements() << " p: " << *this->data_.pressure() << "]";

    return result.str();
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
