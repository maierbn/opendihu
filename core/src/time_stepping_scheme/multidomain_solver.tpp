#include "time_stepping_scheme/multidomain_solver.h"

#include <Python.h>  // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "data_management/multidomain.h"

namespace TimeSteppingScheme
{

template<typename FiniteElementMethodPotentialFlow,typename CellMLAdapterType,typename FiniteElementMethodDiffusion>
MultidomainSolver<FiniteElementMethodPotentialFlow,CellMLAdapterType,FiniteElementMethodDiffusion>::
MultidomainSolver(DihuContext context) :
  TimeSteppingScheme(context["MultidomainSolver"]),
  dataMultidomain_(this->context_), finiteElementMethodPotentialFlow_(this->context_["PotentialFlow"]),
  finiteElementMethodDiffusion_(this->context_["Activation"]), finiteElementMethodDiffusionTotal_(this->context_["Activation"]),
  rankSubset_(std::make_shared<Partition::RankSubset>())
{
  // get python config
  this->specificSettings_ = this->context_.getPythonConfig();

  // parse number of motor units
  nCompartments_ = PythonUtility::getOptionInt(this->specificSettings_, "nCompartments", 1, PythonUtility::NonNegative);
}

template<typename FiniteElementMethodPotentialFlow,typename CellMLAdapterType,typename FiniteElementMethodDiffusion>
void MultidomainSolver<FiniteElementMethodPotentialFlow,CellMLAdapterType,FiniteElementMethodDiffusion>::
advanceTimeSpan()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;

  LOG(DEBUG) << "MultidomainSolver::advanceTimeSpan, timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
    << " n steps: " << this->numberTimeSteps_;

  PetscErrorCode ierr;

  // loop over time steps
  double currentTime = this->startTime_;

  // initialize transmembranePotentialNextTimeStep
  for (int k = 0; k < nCompartments_; k++)
  {
    Vec transmembranePotential = dataMultidomain_.transmembranePotential(k)->getContiguousValuesGlobal();
    Vec transmembraneIncrement = dataMultidomain_.transmembraneIncrementNextTimeStep(k)->getContiguousValuesGlobal();
    cellMLAdapters_[k].evaluateTimesteppingRightHandSideExplicit(transmembranePotential, transmembraneIncrement, 0, currentTime);
  }

  for(int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    if (timeStepNo % this->timeStepOutputInterval_ == 0 && timeStepNo > 0)
    {
      LOG(INFO) << "Multidomain solver, timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }

    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;

    std::vector<Vec> subvectorsRightHandSide(nCompartments_+1);

    // advance cellml
    for (int k = 0; k < nCompartments_; k++)
    {
      // get increment value which was already computed in the previous iteration
      Vec transmembranePotential = dataMultidomain_.transmembranePotential(k)->getContiguousValuesGlobal();
      Vec transmembranePotentialNextTimeStep = dataMultidomain_.transmembranePotentialNextTimeStep(k)->getContiguousValuesGlobal();
      Vec transmembraneIncrement = dataMultidomain_.transmembraneIncrementNextTimeStep(k)->getContiguousValuesGlobal();
      ierr = VecWAXPY(transmembranePotentialNextTimeStep, this->timeStepWidth_, transmembraneIncrement, transmembranePotential); CHKERRV(ierr);
      //  VecWAXPY(Vec w,PetscScalar alpha,Vec x,Vec y), w = alpha x + y

      // compute next delta_u = f(u)
      cellMLAdapters_[k].evaluateTimesteppingRightHandSideExplicit(transmembranePotentialNextTimeStep, transmembraneIncrement, timeStepNo, currentTime);

      // extract ionicCurrentNextTimeStep (-1/Cm I_ion(Vm^(i+1))) from all rates (transmembraneIncrementNextTimeStep)
      std::shared_ptr<FieldVariable::FieldVariable<FunctionSpace,1>> ionicCurrentNextTimeStep
        = dataMultidomain_.ionicCurrentNextTimeStep(k);
      dataMultidomain_.transmembraneIncrementNextTimeStep(k)->extractComponent(0, ionicCurrentNextTimeStep);

      // compute the right hand side entry as rhs[k] = Vm_k^{i} - dt/Cm_k I_ion(Vm_k^{i+1})
      ierr = VecCopy(ionicCurrentNextTimeStep->valuesGlobal(), subvectorsRightHandSide[k]); CHKERRV(ierr);
      ierr = VecScale(subvectorsRightHandSide[k], this->timeStepWidth_); CHKERRV(ierr);
      ierr = VecAXPY(subvectorsRightHandSide[k], 1.0, transmembranePotential); CHKERRV(ierr);
    }

    // fill nested right hand side vector with subvectors
    std::vector<PetscInt> indices(nCompartments_+1);
    std::iota(indices.begin(), indices.end(), 0);
    ierr = VecZeroEntries(subvectorsRightHandSide[nCompartments_]); CHKERRV(ierr);   // set last vector to 0
    ierr = VecNestSetSubVecs(rightHandSide_, nCompartments_+1, indices.data(), subvectorsRightHandSide.data()); CHKERRV(ierr);

    // solve A*u^{t+1} = u^{t} for u^{t+1} where A is the system matrix, solveLinearSystem(b,x)
    this->solveLinearSystem();

    // write the solution from the nested vector back to data
    ierr = VecNestSetSubVecs(rightHandSide_, nCompartments_+1, indices.data(), subvectorsRightHandSide.data()); CHKERRV(ierr);
    for (int k = 0; k < nCompartments_; k++)
    {
      Vec transmembranePotential = dataMultidomain_.transmembranePotential(k)->getContiguousValuesGlobal();
      ierr = VecCopy(subvectorsRightHandSide[k], transmembranePotential); CHKERRV(ierr);
    }
    // get phi_e
    ierr = VecCopy(subvectorsRightHandSide[nCompartments_], dataMultidomain_.extraCellularPotential()->valuesGlobal()); CHKERRV(ierr);

    // stop duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::stop(this->durationLogKey_);

    // write current output values
    this->outputWriterManager_.writeOutput(this->dataMultidomain_, timeStepNo, currentTime);

    // start duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::start(this->durationLogKey_);
  }

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);
}

template<typename FiniteElementMethodPotentialFlow,typename CellMLAdapter,typename FiniteElementMethodDiffusion>
void MultidomainSolver<FiniteElementMethodPotentialFlow,CellMLAdapter,FiniteElementMethodDiffusion>::
run()
{
  // initialize everything
  initialize();

  LOG(DEBUG) << "run potential flow simulation";

  // solve potential flow Laplace problem
  finiteElementMethodPotentialFlow_.run();


  LOG(DEBUG) << "compute gradient field";

  // compute a gradient field from the solution of the potential flow
  dataMultidomain_.flowPotential()->setValues(*finiteElementMethodPotentialFlow_.data().solution());
  dataMultidomain_.flowPotential()->computeGradientField(dataMultidomain_.fibreDirection());

  this->advanceTimeSpan();
}

template<typename FiniteElementMethodPotentialFlow,typename CellMLAdapterType,typename FiniteElementMethodDiffusion>
void MultidomainSolver<FiniteElementMethodPotentialFlow,CellMLAdapterType,FiniteElementMethodDiffusion>::
initialize()
{
  if (this->initialized_)
    return;

  TimeSteppingScheme::initialize();

  LOG(DEBUG) << "initialize multidomain_solver, " << nCompartments_ << " compartments";
  assert(this->specificSettings_);

  // initialize the potential flow finite element method, this also creates the function space
  finiteElementMethodPotentialFlow_.initialize();

  // initialize the data object
  dataMultidomain_.setFunctionSpace(finiteElementMethodPotentialFlow_.functionSpace());
  dataMultidomain_.initialize(nCompartments_);

  // initialize the finite element class, from which only the stiffness matrix is needed
  finiteElementMethodDiffusion_.initialize(dataMultidomain_.fibreDirection());
  finiteElementMethodDiffusionTotal_.initialize(dataMultidomain_.fibreDirection(), nCompartments_);

  initializeCellMLAdapters();

  // parse parameters
  PythonUtility::getOptionVector(this->specificSettings_, "am", nCompartments_, am_);
  PythonUtility::getOptionVector(this->specificSettings_, "cm", nCompartments_, cm_);

  LOG(DEBUG) << "Am: " << am_ << ", Cm: " << cm_;

  // initialize system matrix
  setSystemMatrix(this->timeStepWidth_);

  // initialize linear solver
  if (linearSolver_ == nullptr)
  {
    // retrieve linear solver
    linearSolver_ = this->context_.solverManager()->template solver<Solver::Linear>(
      this->specificSettings_, this->rankSubset_->mpiCommunicator());
  }

  // set matrix used for linear system and preconditioner to ksp context
  assert(this->linearSolver_->ksp());
  PetscErrorCode ierr;
  ierr = KSPSetOperators(*this->linearSolver_->ksp(), systemMatrix_, systemMatrix_); CHKERRV(ierr);

  // initialize rhs vector
  ierr = VecCreate(MPI_COMM_WORLD, &rightHandSide_); CHKERRV(ierr);
  ierr = VecSetType(rightHandSide_, VECNEST); CHKERRV(ierr);

  // initialize solution vector
  ierr = VecDuplicate(rightHandSide_, &solution_); CHKERRV(ierr);

  this->initialized_ = true;
}

template<typename FiniteElementMethodPotentialFlow,typename CellMLAdapterType,typename FiniteElementMethodDiffusion>
void MultidomainSolver<FiniteElementMethodPotentialFlow,CellMLAdapterType,FiniteElementMethodDiffusion>::
initializeCellMLAdapters()
{
  if (VLOG_IS_ON(1))
  {
    VLOG(1) << "CellMLAdapters settings: ";
    PythonUtility::printDict(this->specificSettings_);
  }

  std::vector<PyObject *> cellMLConfigs;
  PythonUtility::getOptionVector(this->specificSettings_, "CellMLAdapters", cellMLConfigs);

  // initialize cellml adapters
  cellMLAdapters_.reserve(nCompartments_);
  if (cellMLConfigs.size() < nCompartments_)
  {
    LOG(FATAL) << "Number of CellMLAdapters (" << cellMLConfigs.size() << ") is smaller than number of compartments (" << nCompartments_ << ")";
  }

  for (int k = 0; k < nCompartments_; k++)
  {

    cellMLAdapters_.emplace_back(this->context_.createSubContext(cellMLConfigs[k]));
    cellMLAdapters_[k].initialize();

    // initialize cellml states
    cellMLAdapters_[k].setInitialValues(dataMultidomain_.transmembranePotential(k));
  }

}

template<typename FiniteElementMethodPotentialFlow,typename CellMLAdapterType,typename FiniteElementMethodDiffusion>
void MultidomainSolver<FiniteElementMethodPotentialFlow,CellMLAdapterType,FiniteElementMethodDiffusion>::
setSystemMatrix(double timeStepWidth)
{
  std::vector<Mat> submatrices(MathUtility::sqr(nCompartments_+1),NULL);

  // fill submatrices, empty submatrices stay be NULL

  PetscErrorCode ierr;
  Mat stiffnessMatrix = finiteElementMethodDiffusion_.data().stiffnessMatrix()->valuesGlobal();

  // set diagonals
  for (int k = 0; k < nCompartments_; k++)
  {
    // right column matrix
    double prefactor = this->timeStepWidth_ / (am_[k]*cm_[k]);
    // create matrix as copy of stiffnessMatrix
    Mat matrixOnRightColumn;
    ierr = MatConvert(stiffnessMatrix, MATSAME, MAT_INITIAL_MATRIX, &matrixOnRightColumn); CHKERRV(ierr);

    // scale with prefactor
    ierr = MatScale(matrixOnRightColumn, prefactor); CHKERRV(ierr);

    // set on right column of the system matrix
    submatrices[k*(nCompartments_+1) + (nCompartments_+1) - 1] = matrixOnRightColumn;

    // ---
    // bottom row matrices
    double f = 1./nCompartments_;
    prefactor = -f;
    // create matrix as copy of stiffnessMatrix
    Mat matrixOnBottomRow;
    ierr = MatConvert(stiffnessMatrix, MATSAME, MAT_INITIAL_MATRIX, &matrixOnBottomRow); CHKERRV(ierr);

    // scale with prefactor
    ierr = MatScale(matrixOnBottomRow, prefactor); CHKERRV(ierr);

    // set on bottom row of the system matrix
    submatrices[((nCompartments_+1) - 1)*(nCompartments_+1) + k] = matrixOnBottomRow;

    // ---
    // diagonal matrix
    // add identity
    Mat matrixOnDiagonalBlock;
    ierr = MatConvert(matrixOnRightColumn, MATSAME, MAT_INITIAL_MATRIX, &matrixOnDiagonalBlock); CHKERRV(ierr);

    ierr = MatShift(matrixOnDiagonalBlock, 1); CHKERRV(ierr);

    // set on diagonal
    submatrices[k*(nCompartments_+1) + k] = matrixOnDiagonalBlock;

    VLOG(2) << "submatrix " << k << ":" << PetscUtility::getStringMatrix(matrixOnDiagonalBlock);
  }

  // set bottom right matrix
  Mat stiffnessMatrixBottomRight = finiteElementMethodDiffusionTotal_.data().stiffnessMatrix()->valuesGlobal();
  ierr = MatScale(stiffnessMatrixBottomRight, -1); CHKERRV(ierr);

  // set on bottom right
  submatrices[(nCompartments_+1)*(nCompartments_+1)-1] = stiffnessMatrixBottomRight;

  // create nested matrix
  ierr = MatCreateNest(this->dataMultidomain_.functionSpace()->meshPartition()->mpiCommunicator(),
                       nCompartments_+1, NULL, nCompartments_+1, NULL, submatrices.data(), &this->systemMatrix_); CHKERRV(ierr);
}

template<typename FiniteElementMethodPotentialFlow,typename CellMLAdapterType,typename FiniteElementMethodDiffusion>
void MultidomainSolver<FiniteElementMethodPotentialFlow,CellMLAdapterType,FiniteElementMethodDiffusion>::
solveLinearSystem()
{
  PetscErrorCode ierr;

  // solve the system, KSPSolve(ksp,b,x)
  ierr = KSPSolve(*this->linearSolver_->ksp(), rightHandSide_, solution_); CHKERRV(ierr);

  int numberOfIterations = 0;
  PetscReal residualNorm = 0.0;
  ierr = KSPGetIterationNumber(*this->linearSolver_->ksp(), &numberOfIterations); CHKERRV(ierr);
  ierr = KSPGetResidualNorm(*this->linearSolver_->ksp(), &residualNorm); CHKERRV(ierr);

  KSPConvergedReason convergedReason;
  ierr = KSPGetConvergedReason(*this->linearSolver_->ksp(), &convergedReason); CHKERRV(ierr);

  VLOG(1) << "Linear system of multidomain problem solved in " << numberOfIterations << " iterations, residual norm " << residualNorm
    << ": " << PetscUtility::getStringLinearConvergedReason(convergedReason);
}

//! return whether the underlying discretizableInTime object has a specified mesh type and is not independent of the mesh type
template<typename FiniteElementMethodPotentialFlow,typename CellMLAdapterType,typename FiniteElementMethodDiffusion>
bool MultidomainSolver<FiniteElementMethodPotentialFlow,CellMLAdapterType,FiniteElementMethodDiffusion>::
knowsMeshType()
{
  return true;
}

} // namespace TimeSteppingScheme
