#include "specialized_solver/solid_mechanics/hyperelasticity/02_hyperelasticity_solver.h"

#include <Python.h>  // has to be the first included header

#include "specialized_solver/solid_mechanics/hyperelasticity/02_petsc_callbacks.h"
#include "solver/nonlinear.h"
#include "control/dihu_context.h"
#include "partition/partitioned_petsc_vec/02_partitioned_petsc_vec_for_hyperelasticity.h"

namespace SpatialDiscretization
{

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticitySolver<Term,withLargeOutput,MeshType,nDisplacementComponents>::
nonlinearSolve()
{
  LOG(TRACE) << "nonlinear solve";
  LOG(DEBUG) << "initial solution: " << combinedVecSolution_->getString();
  // solve the system ∂W_int - ∂W_ext = 0 and J = 1 for displacements and pressure, result will be in solverVariableSolution_, combinedVecSolution_

#ifndef NDEBUG
  this->materialComputeResidual(1.0);   // compute residual with load factor 1.0
  LOG(DEBUG) << "initial residual: " << combinedVecResidual_->getString();
#endif

  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_+std::string("_durationSolve"));

  assert(this->nonlinearSolver_);
  std::shared_ptr<SNES> snes = this->nonlinearSolver_->snes();
  std::shared_ptr<KSP> ksp = this->nonlinearSolver_->ksp();
  this->bestResidualNorm_ = std::numeric_limits<double>::max();
  currentLoadFactor_ = -2*this->loadFactorGiveUpThreshold_;

  // print newline
  LOG(INFO);

  // loop over load loadFactors, i.e. load increments
  std::vector<double> loadFactors(this->loadFactors_);
  for (int loadFactorIndex = 0; loadFactorIndex < loadFactors.size(); loadFactorIndex++)
  {
    double loadFactor = loadFactors[loadFactorIndex];

    previousLoadFactor_ = currentLoadFactor_;
    currentLoadFactor_ = loadFactor;
    if (loadFactors.size() > 1)
    {
      LOG(INFO) << "Nonlinear Solver: load factor " << loadFactor << " of list " << loadFactors;
    }
    if (currentLoadFactor_ < loadFactorGiveUpThreshold_)
    {
      LOG(WARNING) << "Nonlinear solver reached load factor " << currentLoadFactor_ << ", (no. " << loadFactorIndex << ") which "
        << "is below give-up threshold of " << loadFactorGiveUpThreshold_ << ". "
        << "Now abort, use best found solution with residual norm " << bestResidualNorm_;

      // restore best found solution so far
      PetscErrorCode ierr;
      ierr = VecCopy(bestSolution_, solverVariableSolution_); CHKERRV(ierr);
      break;
    }
    else if (currentLoadFactor_ - previousLoadFactor_ < 0 && previousLoadFactor_ - currentLoadFactor_ < loadFactorGiveUpThreshold_)
    {
      LOG(WARNING) << "Nonlinear solver reduced load factor from " << previousLoadFactor_ << " to " << currentLoadFactor_
        << ", (no. " << loadFactorIndex << "), change " << previousLoadFactor_ - currentLoadFactor_ << " "
        << "is below give-up threshold of " << loadFactorGiveUpThreshold_ << ". "
        << "Now abort, use best found solution with residual norm " << bestResidualNorm_;

      // restore best found solution so far
      PetscErrorCode ierr;
      ierr = VecCopy(bestSolution_, solverVariableSolution_); CHKERRV(ierr);
      break;
    }

    // try as many times to solve the nonlinear problem as given in the option nNonlinearSolveCalls
    for (int i = 0; i < this->nNonlinearSolveCalls_; i++)
    {
      LOG(DEBUG) << "------------------  start solve " << i << "/" << this->nNonlinearSolveCalls_ << " ------------------";

      // save initial solution value
      PetscErrorCode ierr;
      if (this->lastSolveSucceeded_)
      {
        ierr = VecCopy(solverVariableSolution_, lastSolution_); CHKERRV(ierr);
      }

      // reset indicator whether the last solve did not encounter a negative jacobian
      this->lastSolveSucceeded_ = true;
      // solve the system nonlinearFunction(displacements) = 0
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

      // compute mean norm
      double normSum = 0;
      for (double norm : this->norms_)
      {
        normSum += norm;
      }
      double meanNorm = normSum / this->norms_.size();
      this->norms_.clear();

      // if the last solution failed, either diverged or got a negative jacobian
      if (!this->lastSolveSucceeded_ || (convergedReason < 0 && residualNorm > meanNorm))
      {
        // add an intermediate load factor
        double lastSuccessfulLoadFactor = 0;
        if (loadFactorIndex > 0)
          lastSuccessfulLoadFactor = loadFactors[loadFactorIndex-1];

        // compute the new load factor to give half the increment size
        double intermediateLoadFactor = 0.5*(currentLoadFactor_ + lastSuccessfulLoadFactor);

        // add new load factor
        loadFactors.insert(loadFactors.begin()+loadFactorIndex, intermediateLoadFactor);
        loadFactorIndex--;

        LOG(INFO) << "Solution failed after " << numberOfIterations << " iterations, residual norm " << residualNorm
          << ", retry with load factor " << intermediateLoadFactor << ": " << PetscUtility::getStringNonlinearConvergedReason(convergedReason) << ", "
          << PetscUtility::getStringLinearConvergedReason(kspConvergedReason);

        this->lastSolveSucceeded_ = false;

        // restore last solution
        ierr = VecCopy(lastSolution_, solverVariableSolution_); CHKERRV(ierr);
      }
      else if (convergedReason < 0 && residualNorm <= meanNorm)
      {
        // if the scheme diverged, but the last residual norm is better than average, accept result as solution
        LOG(INFO) << "Solution accepted in " << numberOfIterations << " iterations, residual norm " << residualNorm
          << " < mean (" << meanNorm << "): " << PetscUtility::getStringNonlinearConvergedReason(convergedReason) << ", "
          << PetscUtility::getStringLinearConvergedReason(kspConvergedReason);
        convergedReason = SNES_CONVERGED_ITS;
      }
      else
      {
        LOG(INFO) << "Solution done in " << numberOfIterations << " iterations, residual norm " << residualNorm
          << ": " << PetscUtility::getStringNonlinearConvergedReason(convergedReason) << ", "
          << PetscUtility::getStringLinearConvergedReason(kspConvergedReason);
      }

      // if the nonlinear scheme converged, finish loop
      if (convergedReason >= 0)
      {
        break;
      }
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

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticitySolver<Term,withLargeOutput,MeshType,nDisplacementComponents>::
postprocessSolution()
{
  // close log file
  this->residualNormLogFile_.close();

  // copy the solution values back to this->data_.displacements() and this->data.pressure() (and this->data_.velocities() for the dynamic case)
  this->setUVP(this->combinedVecSolution_->valuesGlobal());

  // if there are triangles at the corners of the mesh, interpolate the non-dof values
  this->displacementsFunctionSpace_->interpolateNonDofValuesInFieldVariable(this->data_.displacements());
  this->pressureFunctionSpace_->interpolateNonDofValuesInFieldVariable(this->data_.pressure());

  if (nDisplacementComponents == 6)
    this->displacementsFunctionSpace_->interpolateNonDofValuesInFieldVariable(this->data_.velocities());

  // compute the PK2 stress at every node
  this->computePK2StressField();

  LOG(DEBUG) << "solution: " << combinedVecSolution_->getString();

  // update the geometry field by the new displacements, also update field variables for the pressure output writer if needed
  bool usePressureOutputWriter = this->outputWriterManagerPressure_.hasOutputWriters();
  this->data_.updateGeometry(this->displacementsScalingFactor_, usePressureOutputWriter);

  // dump files containing rhs and system matrix
  this->nonlinearSolver_->dumpMatrixRightHandSideSolution(this->solverVariableResidual_, this->solverVariableSolution_);


#ifndef NDEBUG
  checkSolution(solverVariableSolution_);
#endif
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticitySolver<Term,withLargeOutput,MeshType,nDisplacementComponents>::
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
  this->norms_.push_back(currentNorm);

  //T* object = static_cast<T*>(mctx);
  std::stringstream message;
  message << "  Nonlinear solver: iteration " << std::setw(2) << its << ", residual norm " << std::setw(11) << currentNorm;
  if (fabs(experimentalOrderOfConvergence) > 1e-12)
    message << ", e_new=e_old^c with c=" <<  std::setw(3) << experimentalOrderOfConvergence << std::setprecision(6);
  LOG(INFO) << message.str();

  // if we got a better solution and this is load factor 1, store the solution
  if ((currentLoadFactor_ == 1 && currentNorm < bestResidualNorm_) || bestResidualNorm_ == std::numeric_limits<double>::max())
  {
    bestResidualNorm_ = currentNorm;
    PetscErrorCode ierr;
    ierr = VecCopy(solverVariableSolution_, bestSolution_); CHKERRV(ierr);
  }

  static int evaluationNo = 0;  // counter how often this function was called

  if (this->dumpDenseMatlabVariables_)
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
  if (this->residualNormLogFile_.is_open())
  {
    this->residualNormLogFile_ << its << ";" << currentNorm << std::endl;
  }
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticitySolver<Term,withLargeOutput,MeshType,nDisplacementComponents>::
initializePetscCallbackFunctions()
{
  assert(this->nonlinearSolver_);
  std::shared_ptr<SNES> snes = this->nonlinearSolver_->snes();
  std::shared_ptr<KSP> ksp = this->nonlinearSolver_->ksp();

  assert(snes != nullptr);

  // set callback functions that are defined in petsc_callbacks.h
  typedef HyperelasticitySolver<Term,withLargeOutput,MeshType,nDisplacementComponents> ThisClass;

  PetscErrorCode (*callbackNonlinearFunction)(SNES, Vec, Vec, void *)              = *nonlinearFunction<ThisClass>;
  PetscErrorCode (*callbackJacobianAnalytic)(SNES, Vec, Mat, Mat, void *)          = *jacobianFunctionAnalytic<ThisClass>;
  PetscErrorCode (*callbackJacobianFiniteDifferences)(SNES, Vec, Mat, Mat, void *) = *jacobianFunctionFiniteDifferences<ThisClass>;
  PetscErrorCode (*callbackJacobianCombined)(SNES, Vec, Mat, Mat, void *)          = *jacobianFunctionCombined<ThisClass>;
  PetscErrorCode (*callbackMonitorFunction)(SNES, PetscInt, PetscReal, void *)     = *monitorFunction<ThisClass>;

  // set function
  PetscErrorCode ierr;
  ierr = SNESSetFunction(*snes, solverVariableResidual_, callbackNonlinearFunction, this); CHKERRV(ierr);

  // set jacobian
  if (this->useAnalyticJacobian_)
  {
    if (this->useNumericJacobian_)   // use combination of analytic jacobian also with finite differences
    {
      // use the analytic jacobian for the preconditioner and the numeric jacobian (from finite differences) as normal jacobian
      ierr = SNESSetJacobian(*snes, this->solverMatrixAdditionalNumericJacobian_, this->solverMatrixJacobian_, callbackJacobianCombined, this); CHKERRV(ierr);
      //ierr = SNESSetJacobian(*snes, solverMatrixAdditionalNumericJacobian_, solverMatrixAdditionalNumericJacobian_, callbackJacobianCombined, this); CHKERRV(ierr);
      //ierr = SNESSetJacobian(*snes, solverMatrixJacobian_, solverMatrixJacobian_, callbackJacobianCombined, this); CHKERRV(ierr);
      LOG(DEBUG) << "Use combination of numeric and analytic jacobian: " << this->solverMatrixJacobian_;
    }
    else    // use pure analytic jacobian, without fd
    {
      ierr = SNESSetJacobian(*snes, this->solverMatrixJacobian_, this->solverMatrixJacobian_, callbackJacobianAnalytic, this); CHKERRV(ierr);
      LOG(DEBUG) << "Use only analytic jacobian: " << this->solverMatrixJacobian_;
    }
  }
  else
  {
    // set function to compute jacobian from finite differences
    ierr = SNESSetJacobian(*snes, this->solverMatrixJacobian_, this->solverMatrixJacobian_, callbackJacobianFiniteDifferences, this); CHKERRV(ierr);
    LOG(DEBUG) << "Use Finite-Differences approximation for jacobian";
  }

  // prepare log file
  if (this->specificSettings_.hasKey("residualNormLogFilename"))
  {
    std::string logFileName = this->specificSettings_.getOptionString("residualNormLogFilename", "residual_norm.txt");

    this->residualNormLogFile_ = std::ofstream(logFileName, std::ios::out | std::ios::binary | std::ios::trunc);

    if (!this->residualNormLogFile_.is_open())
    {
      LOG(WARNING) << "Could not open log file for residual norm, \"" << logFileName << "\".";
    }
  }

  // set monitor function
  ierr = SNESMonitorSet(*snes, callbackMonitorFunction, this, NULL); CHKERRV(ierr);

}

#if 0
template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticitySolver<Term,withLargeOutput,MeshType,nDisplacementComponents>::
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
#endif

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
bool HyperelasticitySolver<Term,withLargeOutput,MeshType,nDisplacementComponents>::
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
  this->setUVP(solverVariableSolution_);

  // compute the actual output of the nonlinear function
  bool successful = this->materialComputeResidual(currentLoadFactor_);

  // prepare output when there are triangular prism elements at the corners, interpolate values in triangles for dofs that are no real dofs
  combinedVecResidual_->interpolateNonDofValues(this->displacementsFunctionSpace_, this->pressureFunctionSpace_);

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

  return successful;
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
bool HyperelasticitySolver<Term,withLargeOutput,MeshType,nDisplacementComponents>::
evaluateAnalyticJacobian(Vec x, Mat jac)
{
  // copy the values of x to the internal data vectors in this->data_
  this->setUVP(x);

  // compute the jacobian
  return this->materialComputeJacobian();
}

template<typename Term,bool withLargeOutput,typename MeshType,int nDisplacementComponents>
void HyperelasticitySolver<Term,withLargeOutput,MeshType,nDisplacementComponents>::
checkSolution(Vec x)
{
  // Check, if the nonlinear function is zero, i.e., the nonlinear solver has finished with a correct solution.
  // This is only done in debug target

  // evaluate nonlinear function
  evaluateNonlinearFunction(x, solverVariableResidual_);

  // loop over components
  int nComponents = (Term::isIncompressible? 4 : 3);
  for (int componentNo = 0; componentNo < nComponents; componentNo++)
  {
    int nEntries = this->displacementsFunctionSpace_->nDofsLocalWithoutGhosts();
    std::vector<double> values;

    // get all raw values of the component
    if (componentNo < 3)
    {
      values.resize(nEntries);
      combinedVecResidual_->getValues(componentNo, nEntries, this->displacementsFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data());
    }
    else if (componentNo == 3 && Term::isIncompressible)
    {
      nEntries = this->pressureFunctionSpace_->nDofsLocalWithoutGhosts();
      values.resize(nEntries);
      combinedVecResidual_->getValues(componentNo, nEntries, this->pressureFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data());
    }

    // check if values are zero
    for (int i = 0; i < nEntries; i++)
    {
      double value = values[i];
      if (fabs(value) > 1e-5)
        LOG(DEBUG) << "Component " << componentNo << " entry " << i << ", residual is " << value << ", should be zero.";
    }
  }

  solverVariableResidual_ = combinedVecResidual_->valuesGlobal();

  // compute L2 norm
  PetscReal l2NormResidual;
  VecNorm(solverVariableResidual_, NORM_2, &l2NormResidual);
  LOG(DEBUG) << "L2-norm residual: " << l2NormResidual;

  // print result
  if (l2NormResidual < 1e-5)
  {
    LOG(DEBUG) << ANSI_COLOR_GREEN "Root found." ANSI_COLOR_RESET;
  }
  else
  {
    LOG(DEBUG) << ANSI_COLOR_RED "Root NOT found." ANSI_COLOR_RESET;
  }
}

} // namespace TimeSteppingScheme
