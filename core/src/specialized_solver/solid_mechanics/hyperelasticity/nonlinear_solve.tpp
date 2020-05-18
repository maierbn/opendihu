#include "specialized_solver/solid_mechanics/hyperelasticity/hyperelasticity_solver.h"

#include <Python.h>  // has to be the first included header

#include "specialized_solver/solid_mechanics/hyperelasticity/petsc_callbacks.h"
#include "solver/nonlinear.h"
#include "control/dihu_context.h"
#include "partition/partitioned_petsc_vec/02_partitioned_petsc_vec_for_hyperelasticity.h"

namespace SpatialDiscretization
{

template<typename Term,typename MeshType, int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
nonlinearSolve()
{
  LOG(TRACE) << "nonlinear solve";
  LOG(DEBUG) << "initial solution: " << combinedVecSolution_->getString();
  // solve the system ∂W_int - ∂W_ext = 0 and J = 1 for displacements and pressure, result will be in solverVariableSolution_, combinedVecSolution_

#ifndef NDEBUG
  materialComputeResidual(1.0);   // compute residual with load factor 1.0
  LOG(DEBUG) << "initial residual: " << combinedVecResidual_->getString();
#endif

  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_+std::string("_durationSolve"));

  assert(nonlinearSolver_);
  std::shared_ptr<SNES> snes = nonlinearSolver_->snes();
  std::shared_ptr<KSP> ksp = nonlinearSolver_->ksp();

  // newline
  LOG(INFO);

  // loop over load loadFactors, i.e. load increments
  for (double loadFactor: loadFactors_)
  {
    currentLoadFactor_ = loadFactor;
    if (loadFactors_.size() > 1)
    {
      LOG(INFO) << "Nonlinear Solver: load factor " << loadFactor << " of list " << loadFactors_;
    }


    // try two times to solve the nonlinear problem
    for (int i = 0; i < nNonlinearSolveCalls_; i++)
    {
      LOG(DEBUG) << "------------------  start solve " << i << "/" << nNonlinearSolveCalls_ << " ------------------";

      // solve the system nonlinearFunction(displacements) = 0
      PetscErrorCode ierr;
      ierr = SNESSolve(*snes, NULL, solverVariableSolution_); CHKERRV(ierr);

      // get information about the solution process
      PetscInt numberOfIterations = 0;
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

      // if the nonlinear scheme converged, finish loop
      if (convergedReason >= 0)
        break;
    }

    // reset value of last residual norm that is needed for computational of experimental order of convergence
    lastNorm_ = 0;
    secondLastNorm_ = 0;

    // write current output values
    if (this->outputWriterManagerLoadIncrements_.hasOutputWriters())
    {
      this->outputWriterManagerLoadIncrements_.writeOutput(this->data_, 1, endTime_);
    }
  }

  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_+std::string("_durationSolve"));

}

template<typename Term,typename MeshType, int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
postprocessSolution()
{
  // close log file
  residualNormLogFile_.close();

  // copy the solution values back to this->data_.displacements() and this->data.pressure() (and this->data_.velocities() for the dynamic case)
  setUVP(combinedVecSolution_->valuesGlobal());

  // compute the PK2 stress at every node
  computePK2StressField();

  LOG(DEBUG) << "solution: " << combinedVecSolution_->getString();

  // update the geometry field by the new displacements, also update field variables for the pressure output writer if needed
  bool usePressureOutputWriter = this->outputWriterManagerPressure_.hasOutputWriters();
  this->data_.updateGeometry(displacementsScalingFactor_, usePressureOutputWriter);

  // dump files containing rhs and system matrix
  nonlinearSolver_->dumpMatrixRightHandSideSolution(solverVariableResidual_, solverVariableSolution_);


#ifndef NDEBUG
  checkSolution(solverVariableSolution_);
#endif
}

template<typename Term,typename MeshType, int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
monitorSolvingIteration(SNES snes, PetscInt its, PetscReal currentNorm)
{
  // compute experimental order of convergence which is a measure for the current convergence velocity
  //PetscReal experimentalOrderOfConvergence = 0;

  //if (secondLastNorm_ != 0)
  //  experimentalOrderOfConvergence = log(lastNorm_ / currentNorm) / log(secondLastNorm_ / lastNorm_);

  // derivation of experimentalOrderOfConvergence
  // e_current = e_old ^ c = exp(c*log(e_old)) => c = log(e_current) / log(e_old)
  PetscReal experimentalOrderOfConvergence = log(currentNorm) / log(lastNorm_);

  secondLastNorm_ = lastNorm_;
  lastNorm_ = currentNorm;

  //T* object = static_cast<T*>(mctx);
  std::stringstream message;
  message << "  Nonlinear solver: iteration " << std::setw(2) << its << ", residual norm " << std::setw(11) << currentNorm
          << ", e_new=e_old^c with c=" <<  std::setw(3) << experimentalOrderOfConvergence << std::setprecision(6);
  LOG(INFO) << message.str();

  static int evaluationNo = 0;  // counter how often this function was called

  if (dumpDenseMatlabVariables_)
  {
    // dump input vector
    // dumpVector(std::string filename, std::string format, Vec &vector, MPI_Comm mpiCommunicator, int componentNo=0, int nComponents=1);
    std::stringstream filename;
    filename << "out/x" << std::setw(3) << std::setfill('0') << evaluationNo;
    //PetscUtility::dumpVector(filename.str(), "matlab", solverVariableSolution_, displacementsFunctionSpace->meshPartition()->mpiCommunicator());
    combinedVecSolution_->dumpGlobalNatural(filename.str());

    filename.str("");
    filename << "out/r" << std::setw(3) << std::setfill('0') << evaluationNo;
    //PetscUtility::dumpVector(filename.str(), "matlab", solverVariableSolution_, displacementsFunctionSpace->meshPartition()->mpiCommunicator());
    combinedVecResidual_->dumpGlobalNatural(filename.str());
  }

  evaluationNo++;

  // if log file was given, write residual norm to log file
  if (residualNormLogFile_.is_open())
  {
    residualNormLogFile_ << its << ";" << currentNorm << std::endl;
  }
}

template<typename Term,typename MeshType, int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
debug()
{
  materialComputeInternalVirtualWork();
  // now, solverVariableResidual_, which is the globalValues() of combinedVecResidual_, contains δW_int
  // also the pressure equation residual has been set at the last component

  // for static case: F = δW_int - δW_ext

  // compute F = δW_int - δW_ext,
  // δW_ext = int_∂Ω T_a phi_L dS was precomputed in initialize (materialComputeExternalVirtualWorkDead()), in variable externalVirtualWorkDead_
  // for static case, externalVirtualWorkDead_ = externalVirtualWorkDead_
  PetscErrorCode ierr;
  ierr = VecAXPY(solverVariableResidual_, -1.0, externalVirtualWorkDead_); CHKERRV(ierr);

  LOG(DEBUG) << "check if initial solution matches static problem. Residual:   " << getString(solverVariableResidual_);

  // compute the residual norm

  PetscReal l2NormResidual;
  VecNorm(solverVariableResidual_, NORM_2, &l2NormResidual);
  LOG(DEBUG) << "L2-norm residual: " << l2NormResidual;

  LOG(FATAL) << "done";

  return;


  PetscInt nRows, nColumns;
  MatGetLocalSize(solverMatrixAdditionalNumericJacobian_, &nRows, &nColumns);

  std::vector<double> entriesNumericJacobian;
  PetscUtility::getMatrixEntries(solverMatrixAdditionalNumericJacobian_, entriesNumericJacobian);

  std::vector<double> entriesAnalyticJacobian;
  PetscUtility::getMatrixEntries(solverMatrixJacobian_, entriesAnalyticJacobian);

  LOG(DEBUG) << "nRows: " << nRows << ", nColumns: " << nColumns << ", nEntries: " << entriesNumericJacobian.size() << "=" << entriesAnalyticJacobian.size() << " should be " << nRows*nColumns;

  bool uSame = true;
  double maxError = 0;
  for (int rowNo = 0; rowNo < nRows; rowNo++)
  {
    for (int columnNo = 0; columnNo < nColumns; columnNo++)
    {
      double error = entriesNumericJacobian[rowNo*nColumns + columnNo] - entriesAnalyticJacobian[rowNo*nColumns + columnNo];
      if (fabs(error) > maxError)
        maxError = error;
      if (fabs(error) > 1e-5)
      {
        uSame = false;
        LOG(DEBUG) << "  (" << rowNo << "," << columnNo << ") numeric: " << entriesNumericJacobian[rowNo*nColumns + columnNo]
          << ", analytic: " << entriesAnalyticJacobian[rowNo*nColumns + columnNo];
        //break;
      }
    }
  }
/*
  bool pSame = true;
  for (int rowNo = 81-12; rowNo < 89-12; rowNo++)
  {
    for (int columnNo = 0; columnNo < 89-12; columnNo++)
    {
      if (fabs(entriesNumericJacobian[rowNo*nColumns + columnNo] - entriesAnalyticJacobian[rowNo*nColumns + columnNo]) > 1e-5)
      {
        pSame = false;
        LOG(DEBUG) << " p (" << rowNo << "," << columnNo << ") numeric: " << entriesNumericJacobian[rowNo*nColumns + columnNo]
          << ", analytic: " << entriesAnalyticJacobian[rowNo*nColumns + columnNo];
        //break;
      }
    }
  }*/

  MPI_Barrier(this->combinedVecSolution_->meshPartition()->mpiCommunicator());
  LOG(DEBUG) << "same: " << uSame << ", max error: " << maxError;// << ", pSame: " << pSame;

  return;


  //materialComputeResidual();   // solverVariableSolution_ -> solverVariableResidual_ respective combinedVecSolution_ -> combinedVecResidual_

  PetscInt i = 46;  // 46  17
  PetscInt j = 47;  // 47  18
  double epsilon = 1e-10;

  // compute F(x) - F(x+eps_i)
  Vec *x;
  Vec *f, *diff;
  VecDuplicateVecs(solverVariableSolution_, 5, &x);
  VecDuplicateVecs(solverVariableSolution_, 5, &f);
  VecDuplicateVecs(solverVariableSolution_, 2, &diff);

  PetscInt begin,end;
  VecGetOwnershipRange(solverVariableSolution_,&begin,&end);

  // set x0 = x
  VecCopy(solverVariableSolution_, x[0]);

  // set x1 = x+eps_i
  VecCopy(solverVariableSolution_, x[1]);

  // set x2 = x+eps_j
  VecCopy(solverVariableSolution_, x[2]);

  if (i >= begin && i < end)
  {
    LOG(DEBUG) << "set epsilon at i=" << i;
    VecSetValue(x[1], i, epsilon, ADD_VALUES);
  }
  VecAssemblyBegin(x[1]);
  VecAssemblyEnd(x[1]);

  // compute diff0 = F(x1) - f(x0)
  evaluateNonlinearFunction(x[0], f[0]);
  evaluateNonlinearFunction(x[1], f[1]);

  VecWAXPY(diff[0],-1.0,f[1],f[0]);  // WAXPY(w,a,x,y), w = a*x + y
  VecScale(diff[0], 1./epsilon);

  LOG(DEBUG) << "(F(x+eps_i) - F(x)) / eps = diff0 = " << PetscUtility::getStringVector(diff[0]);


  if (j >= begin && j < end)
  {
    LOG(DEBUG) << "set epsilon at j=" << j;
    VecSetValue(x[2], j, epsilon, ADD_VALUES);
  }
  VecAssemblyBegin(x[2]);
  VecAssemblyEnd(x[2]);

  // compute diff1 = F(x2) - f(x0)
  evaluateNonlinearFunction(x[2], f[2]);

  VecWAXPY(diff[1],-1.0,f[2],f[0]);  // WAXPY(w,a,x,y), w = a*x + y
  VecScale(diff[1], 1./epsilon);

  LOG(DEBUG) << "(F(x+eps_j) - F(x)) / eps = diff1 = " << PetscUtility::getStringVector(diff[1]);

  if (j >= begin && j < end)
  {
    double value;
    VecGetValues(diff[0],1,&j,&value);
    LOG(DEBUG) << "j = " << j << ", dF_j(x)/dx_i = " << value;
  }
  if (i >= begin && i < end)
  {
    double value;
    VecGetValues(diff[1],1,&i,&value);
    LOG(DEBUG) << "i = " << i << ", dF_i(x)/dx_j = " << value;
  }

  MPI_Barrier(combinedVecSolution_->meshPartition()->mpiCommunicator());
  LOG(FATAL) << "end";
}

template<typename Term,typename MeshType, int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
initializeSolutionVariable()
{
  // set variable to all zero and dirichlet boundary condition value
  LOG(DEBUG) << "zeroEntries, representation: " << combinedVecSolution_->currentRepresentation();
  combinedVecSolution_->zeroEntries();

  LOG(DEBUG) << "values: " << PetscUtility::getStringVector(solverVariableSolution_);
  LOG(DEBUG) << "after initialization: " << combinedVecSolution_->getString();
}

template<typename Term,typename MeshType, int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
initializePetscCallbackFunctions()
{
  assert(nonlinearSolver_);
  std::shared_ptr<SNES> snes = nonlinearSolver_->snes();
  std::shared_ptr<KSP> ksp = nonlinearSolver_->ksp();

  assert(snes != nullptr);
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
  typedef HyperelasticitySolver ThisClass;

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
      // use the analytic jacobian for the preconditioner and the numeric jacobian (from finite differences) as normal jacobian
      ierr = SNESSetJacobian(*snes, solverMatrixAdditionalNumericJacobian_, solverMatrixJacobian_, callbackJacobianCombined, this); CHKERRV(ierr);
      //ierr = SNESSetJacobian(*snes, solverMatrixAdditionalNumericJacobian_, solverMatrixAdditionalNumericJacobian_, callbackJacobianCombined, this); CHKERRV(ierr);
      //ierr = SNESSetJacobian(*snes, solverMatrixJacobian_, solverMatrixJacobian_, callbackJacobianCombined, this); CHKERRV(ierr);
      LOG(DEBUG) << "Use combination of numeric and analytic jacobian: " << solverMatrixJacobian_;
    }
    else    // use pure analytic jacobian, without fd
    {
      ierr = SNESSetJacobian(*snes, solverMatrixJacobian_, solverMatrixJacobian_, callbackJacobianAnalytic, this); CHKERRV(ierr);
      LOG(DEBUG) << "Use only analytic jacobian: " << solverMatrixJacobian_;
    }
  }
  else
  {
    // set function to compute jacobian from finite differences
    ierr = SNESSetJacobian(*snes, solverMatrixJacobian_, solverMatrixJacobian_, callbackJacobianFiniteDifferences, this); CHKERRV(ierr);
    LOG(DEBUG) << "Use Finite-Differences approximation for jacobian";
  }

  // prepare log file
  if (this->specificSettings_.hasKey("residualNormLogFilename"))
  {
    std::string logFileName = this->specificSettings_.getOptionString("residualNormLogFilename", "residual_norm.txt");

    residualNormLogFile_ = std::ofstream(logFileName, std::ios::out | std::ios::binary | std::ios::trunc);

    if (!residualNormLogFile_.is_open())
    {
      LOG(WARNING) << "Could not open log file for residual norm, \"" << logFileName << "\".";
    }
  }

  // set monitor function
  ierr = SNESMonitorSet(*snes, callbackMonitorFunction, this, NULL); CHKERRV(ierr);

}

template<typename Term,typename MeshType, int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
evaluateNonlinearFunction(Vec x, Vec f)
{
  //VLOG(1) << "evaluateNonlinearFunction at " << getString(x);

  // determine if Vecs need to be backed up and instead x,f should take the place of solverVariableSolution_,solverVariableResidual_
  // this happens in computation of the numeric jacobian
  bool backupVecs = f != solverVariableResidual_ || x != solverVariableSolution_;

  //VLOG(1) << "backupVecs: " << backupVecs;

  Vec backupSolution;
  Vec backupResidual;
  if (backupVecs)
  {
    // backup the Vecs of solverVariableSolution_ and solverVariableResidual_
    backupSolution = combinedVecSolution_->valuesGlobal();
    backupResidual = combinedVecResidual_->valuesGlobal();
    
    // assign x and f to the variables solverVariableSolution_ and solverVariableResidual_
    combinedVecSolution_->valuesGlobalReference() = x;
    solverVariableSolution_ = x;
    
    combinedVecResidual_->valuesGlobalReference() = f;
    solverVariableResidual_ = f;

    // the following would have done the same, but gives the following error, when Petsc is compiled in debug mode:
    // PETSC ERROR: VecSetErrorIfLocked() line 556 in /store/software/opendihu/dependencies/petsc/src/petsc-3.12.3/include/petscvec.h  Vec is already locked for read-only or read/write access, argument # 1
    //VecSwap(x, solverVariableSolution_);
    //VecSwap(f, solverVariableResidual_);
  }

  // set the solverVariableSolution_ values in displacements, velocities and pressure, this is needed for materialComputeResidual
  setUVP(solverVariableSolution_);

  // compute the actual output of the nonlinear function
  materialComputeResidual(currentLoadFactor_);

  //VLOG(1) << "solverVariableResidual_: " << combinedVecResidual_->getString();

  // restore the values of solverVariableSolution_ and solverVariableResidual_ to their original pointer
  if (backupVecs)
  {
    combinedVecSolution_->valuesGlobalReference() = backupSolution;
    combinedVecResidual_->valuesGlobalReference() = backupResidual;
    solverVariableSolution_ = combinedVecSolution_->valuesGlobal();
    solverVariableResidual_ = combinedVecResidual_->valuesGlobal();

    //VecSwap(x, solverVariableSolution_);
    //VecSwap(f, solverVariableResidual_);
  }
  //VLOG(1) << "f: " << getString(f);
}

template<typename Term,typename MeshType, int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
evaluateAnalyticJacobian(Vec x, Mat jac)
{
  // copy the values of x to the internal data vectors in this->data_
  setUVP(x);

  // compute the jacobian
  materialComputeJacobian();

  //MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY);
  //MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY);
}

template<typename Term,typename MeshType, int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
setDisplacementsAndPressureFromCombinedVec(Vec x, std::shared_ptr<DisplacementsFieldVariableType> u,
                                           std::shared_ptr<PressureFieldVariableType> p)
{
  assert(nDisplacementComponents == 3);

  // copy entries of combined vector x to this->data_.displacements() and this->data_.pressure()
  std::vector<double> values;

  if (VLOG_IS_ON(1))
  {
    PetscUtility::getVectorEntries(x, values);
    VLOG(1) << "setDisplacementsAndPressureFromCombinedVec, x=" << PetscUtility::getStringVector(x);
  }

  bool backupVecs = false;
  if (x != solverVariableSolution_)
  {
    backupVecs = true;
  }

  // move the values of x to variable solverVariableSolution_
  if (backupVecs)
  {
    // this happens in computation of the numeric jacobian
    VecSwap(x, solverVariableSolution_);
  }

  if (!u && !p)
  {
    u = this->data_.displacements();
    if (Term::isIncompressible)
      p = this->data_.pressure();   // p is only needed for incompressible formulation
  }

  // set displacement entries
  u->zeroGhostBuffer();
  if (p)
  {
    p->zeroGhostBuffer();
  }
  for (int componentNo = 0; componentNo < 3; componentNo++)
  {
    int nEntries = displacementsFunctionSpace_->nDofsLocalWithoutGhosts();
    values.resize(nEntries);
    combinedVecSolution_->getValues(componentNo, nEntries, displacementsFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data());

    u->setValues(componentNo, nEntries, displacementsFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data());
  }
  u->finishGhostManipulation();
  u->startGhostManipulation();

  // set pressure entries
  if (p)
  {
    int nEntries = pressureFunctionSpace_->nDofsLocalWithoutGhosts();
    values.resize(nEntries);

    const int pressureComponent = nDisplacementComponents;
    combinedVecSolution_->getValues(pressureComponent, nEntries, pressureFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data());
    p->setValues(0, nEntries, pressureFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data());
    p->finishGhostManipulation();
    p->startGhostManipulation();
  }

  // undo the backup operation
  if (backupVecs)
  {
    VecSwap(x, solverVariableSolution_);
  }

  //VLOG(1) << *u;
  //VLOG(1) << *p;
}

template<typename Term,typename MeshType, int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
setDisplacementsVelocitiesAndPressureFromCombinedVec(Vec x,
                                                     std::shared_ptr<DisplacementsFieldVariableType> u,
                                                     std::shared_ptr<DisplacementsFieldVariableType> v,
                                                     std::shared_ptr<PressureFieldVariableType> p)
{
  assert(nDisplacementComponents == 6);

  // copy entries of combined vector x to this->data_.displacements() and this->data_.pressure()
  std::vector<double> values;

  if (VLOG_IS_ON(1))
  {
    PetscUtility::getVectorEntries(x, values);
    VLOG(1) << "setDisplacementsVelocitiesAndPressureFromCombinedVec, x=" << PetscUtility::getStringVector(x);
  }

  bool backupVecs = false;
  if (x != solverVariableSolution_)
  {
    backupVecs = true;
  }

  // move the values of x to variable solverVariableSolution_
  if (backupVecs)
  {
    // this happens in computation of the numeric jacobian
    VecSwap(x, solverVariableSolution_);
  }

  if (!u && !v && !p)
  {
    u = this->data_.displacements();
    v = this->data_.velocities();

    if (Term::isIncompressible)
      p = this->data_.pressure();   // p is only needed for incompressible formulation
  }

  // set displacement entries
  u->zeroGhostBuffer();
  for (int componentNo = 0; componentNo < 3; componentNo++)
  {
    int nEntries = displacementsFunctionSpace_->nDofsLocalWithoutGhosts();
    values.resize(nEntries);
    combinedVecSolution_->getValues(componentNo, nEntries, displacementsFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data());

    //if (VLOG_IS_ON(1))
      LOG(DEBUG) << "setDisplacementsVelocitiesAndPressureFromCombinedVec, " << nEntries << " u values: " << values;

    u->setValues(componentNo, nEntries, displacementsFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data());
  }
  u->finishGhostManipulation();
  u->startGhostManipulation();

  // set velocity entries
  if (v)
  {
    v->zeroGhostBuffer();
    for (int componentNo = 0; componentNo < 3; componentNo++)
    {
      int nEntries = displacementsFunctionSpace_->nDofsLocalWithoutGhosts();
      values.resize(nEntries);
      combinedVecSolution_->getValues(3 + componentNo, nEntries, displacementsFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data());

    if (VLOG_IS_ON(1))
      VLOG(1) << "setDisplacementsVelocitiesAndPressureFromCombinedVec, " << nEntries << " v values: " << values;

      v->setValues(componentNo, nEntries, displacementsFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data());
    }
    v->finishGhostManipulation();
    v->startGhostManipulation();
  }
  // set pressure entries
  if (p)
  {
    p->zeroGhostBuffer();
    int nEntries = pressureFunctionSpace_->nDofsLocalWithoutGhosts();
    values.resize(nEntries);

    const int pressureComponent = nDisplacementComponents;
    combinedVecSolution_->getValues(pressureComponent, nEntries, pressureFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data());

    if (VLOG_IS_ON(1))
      VLOG(1) << "setDisplacementsVelocitiesAndPressureFromCombinedVec, " << nEntries << " p values: " << values;

    p->setValues(0, nEntries, pressureFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data());
    p->finishGhostManipulation();
    p->startGhostManipulation();
  }

  // undo the backup operation
  if (backupVecs)
  {
    VecSwap(x, solverVariableSolution_);
  }

  //VLOG(1) << *u;
  //VLOG(1) << *p;
}

template<typename Term,typename MeshType, int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
setUVP(Vec x)
{
  if (nDisplacementComponents == 3)
  {
    setDisplacementsAndPressureFromCombinedVec(x);
  }
  else if (nDisplacementComponents == 6)
  {
    setDisplacementsVelocitiesAndPressureFromCombinedVec(x);
  }
}

template<typename Term,typename MeshType, int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
dumpJacobianMatrix(Mat jac)
{
  if (!dumpDenseMatlabVariables_)
    return;

  static int evaluationNo = 0;  // counter how often this function was called

  std::stringstream filename;
  filename << "out/jac_" << std::setw(3) << std::setfill('0') << evaluationNo << "r" <<  DihuContext::nRanksCommWorld();

  evaluationNo++;

  // if there are both numeric and analytic jacobian in use, this is the numeric jacobian
  if (jac == solverMatrixAdditionalNumericJacobian_)
  {
    filename << "_numeric";
    combinedMatrixAdditionalNumericJacobian_->dumpMatrixGlobalNatural(filename.str());
  }
  else if (jac == solverMatrixJacobian_)
  {
    // this is the normal jacobian, either numeric or analytic, if only one of both is in use

    // if both matrices are used
    if (solverMatrixAdditionalNumericJacobian_ != PETSC_NULL)
    {
      filename << "_analytic";
      combinedMatrixJacobian_->dumpMatrixGlobalNatural(filename.str());

      evaluationNo--;

      // if matrix is not all zeros
      double numericJacobianNorm = 0;
      MatNorm(solverMatrixAdditionalNumericJacobian_, NORM_1, &numericJacobianNorm);

      LOG(DEBUG) << "numericJacobianNorm: " << numericJacobianNorm;

      if (numericJacobianNorm > 1e-5)
      {
        // compute difference between analytic and numeric jacobian
        Mat difference;
        PetscErrorCode ierr;
        ierr = MatDuplicate(solverMatrixAdditionalNumericJacobian_, MAT_COPY_VALUES, &difference); CHKERRV(ierr);
        ierr = MatCopy(solverMatrixAdditionalNumericJacobian_, difference, SAME_NONZERO_PATTERN); CHKERRV(ierr);
        MatAXPY(difference, -1, solverMatrixJacobian_, DIFFERENT_NONZERO_PATTERN);

        double norm1 = 0;
        double normF = 0;
        double normInf = 0;
        MatNorm(difference, NORM_1, &norm1);
        MatNorm(difference, NORM_FROBENIUS, &normF);
        MatNorm(difference, NORM_INFINITY, &normInf);
        LOG(INFO) << "difference between analytic and numeric jacobian matrices: "
          << "1-norm: " << norm1 << ", frobenius norm: " << normF << ", infinity norm: " << normInf;

        // compute differences for submatrices
        int nRows = 2;
        if (nDisplacementComponents == 6)
          nRows = 3;
        if (!Term::isIncompressible)    // compressible formulation does not have pressure component
          nRows--;

        for (int i = 0; i < nRows; i++)
        {
          for (int j = 0; j < nRows; j++)
          {
            Mat analyticJacobianSubmatrix = combinedMatrixJacobian_->getSubmatrix(i,j);
            Mat numericJacobianSubmatrix = combinedMatrixAdditionalNumericJacobian_->getSubmatrix(i,j);

            Mat differenceSubmatrix;
            ierr = MatDuplicate(analyticJacobianSubmatrix, MAT_COPY_VALUES, &differenceSubmatrix); CHKERRV(ierr);
            ierr = MatCopy(analyticJacobianSubmatrix, differenceSubmatrix, SAME_NONZERO_PATTERN); CHKERRV(ierr);
            MatAXPY(differenceSubmatrix, -1, numericJacobianSubmatrix, DIFFERENT_NONZERO_PATTERN);

            MatNorm(differenceSubmatrix, NORM_1, &norm1);
            MatNorm(differenceSubmatrix, NORM_FROBENIUS, &normF);
            MatNorm(differenceSubmatrix, NORM_INFINITY, &normInf);
            LOG(INFO) << " submatrix (" << i << "," << j << "): 1-norm: " << norm1 << ", frobenius norm: " << normF << ", infinity norm: " << normInf;
          }
        }

        if (norm1 > 1)
          LOG(ERROR) << "norm mismatch";
      }
    }
    else
    {
      combinedMatrixJacobian_->dumpMatrixGlobalNatural(filename.str());
    }
  }
  else
  {
    LOG(ERROR) << "Could not output jacobian matrix " << jac << " (solverMatrixJacobian_: " << solverMatrixJacobian_ << ", solverMatrixAdditionalNumericJacobian_: " << solverMatrixAdditionalNumericJacobian_ << ")";
  }
}

template<typename Term,typename MeshType, int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
checkSolution(Vec x)
{
  // check if function is zero
  evaluateNonlinearFunction(x, solverVariableResidual_);

  int nComponents = (Term::isIncompressible? 4 : 3);

  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    int nEntries = displacementsFunctionSpace_->nDofsLocalWithoutGhosts();
    std::vector<double> values;

    if (componentNo < 3)
    {
      values.resize(nEntries);
      combinedVecResidual_->getValues(componentNo, nEntries, displacementsFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data());
    }
    else if (componentNo == 3 && Term::isIncompressible)
    {
      nEntries = pressureFunctionSpace_->nDofsLocalWithoutGhosts();
      values.resize(nEntries);
      combinedVecResidual_->getValues(componentNo, nEntries, pressureFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data());
    }

    for (int i = 0; i < nEntries; i++)
    {
      double value = values[i];
      if (fabs(value) > 1e-5)
        LOG(DEBUG) << "Component " << componentNo << " entry " << i << ", residual is " << value << ", should be zero.";
    }
  }

  solverVariableResidual_ = combinedVecResidual_->valuesGlobal();

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

template<typename Term,typename MeshType, int nDisplacementComponents>
std::string HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
getString(Vec x)
{
  if (x == solverVariableSolution_)
  {
    return combinedVecSolution_->getString();
  }
  else if (x == solverVariableResidual_)
  {
    return combinedVecResidual_->getString();
  }
  else if (x == externalVirtualWorkDead_)
  {
    return combinedVecExternalVirtualWorkDead_->getString();
  }
  else
  {
    LOG(FATAL) << "this should not be called";
  }
  return std::string("no getString representation");
}

} // namespace TimeSteppingScheme
