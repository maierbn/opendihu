#include "time_stepping_scheme/multidomain_solver.h"

#include <Python.h>  // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "data_management/multidomain.h"

#define MONODOMAIN

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

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_, this->specificSettings_);

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

  // loop over time steps
  for(int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    if (timeStepNo % this->timeStepOutputInterval_ == 0 && timeStepNo > 0)
    {
      LOG(INFO) << "Multidomain solver, timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime;
    }

    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;

    VLOG(1) << "Godunov splitting";
    VLOG(1) << "---- reaction term";

    // Godunov splitting
    // advance cellml
    for (int k = 0; k < nCompartments_; k++)
    {
      VLOG(1) << "k=" << k;

      // extract transmembranePotential from subcellularStates (TODO: could this be done using globalValues(0)?)
      dataMultidomain_.subcellularStates(k)->extractComponent(0, dataMultidomain_.transmembranePotential(k));

      // get subcellular variable vectors
      Vec subcellularStates = dataMultidomain_.subcellularStates(k)->getContiguousValuesGlobal();
      Vec subcellularIncrement = dataMultidomain_.subcellularIncrement(k)->getContiguousValuesGlobal();

      // compute next delta_y = f(y)
      cellMLAdapters_[k].evaluateTimesteppingRightHandSideExplicit(subcellularStates, subcellularIncrement, timeStepNo, currentTime);

      // compute next subcellular states (y += dt*f(y) ), now subcellularStates[0] contains Vm_k^{(*)} = Vm_k^{i} - dt/Cm_k*I_ion(Vm_k^{i})
      ierr = VecAXPY(subcellularStates, this->timeStepWidth_, subcellularIncrement); CHKERRV(ierr);

      //VLOG(2) << "subcellularStates: " << *dataMultidomain_.subcellularStates(k);
      //VLOG(2) << "subcellularIncrement: " << *dataMultidomain_.subcellularIncrement(k);

      // extract ionicCurrent (-1/Cm I_ion(Vm^(i+1))) from all rates (subcellularIncrement)
      //dataMultidomain_.subcellularIncrement(k)->extractComponent(0, dataMultidomain_.ionicCurrent(k));
      //VLOG(2) << "ionicCurrent (should be first component of subcellularIncrement): " << *ionicCurrent;

      // compute the right hand side entry as rhs[k] = Vm_k^{*} = Vm_k^{i} - dt/Cm_k*I_ion(Vm_k^{i})
      // extract the 0th component of the subcellular states which is Vm at the intermediate timestep
      dataMultidomain_.subcellularStates(k)->extractComponent(0, dataMultidomain_.transmembranePotential(k));
      ierr = VecCopy(dataMultidomain_.transmembranePotential(k)->valuesLocal(), subvectorsRightHandSide_[k]); CHKERRV(ierr);   /// rhs[k] = Vm_k^{*}

      //VLOG(2) << "dt = " << this->timeStepWidth_;
      //VLOG(2) << "k=" << k << ", rhs: " << PetscUtility::getStringVector(subvectorsRightHandSide_[k]);
    }

    LOG(DEBUG) << " Vm: ";
    dataMultidomain_.subcellularStates(0)->extractComponent(0, dataMultidomain_.transmembranePotential(0));
    LOG(DEBUG) << *dataMultidomain_.transmembranePotential(0);

    // advance diffusion
    VLOG(1) << "---- diffusion term";

    // fill nested right hand side vector with subvectors
    std::vector<PetscInt> indices(nCompartments_+1);
    std::iota(indices.begin(), indices.end(), 0);
    ierr = VecNestSetSubVecs(rightHandSide_, nCompartments_+1, indices.data(), subvectorsRightHandSide_.data()); CHKERRV(ierr);

    // solve A*u^{t+1} = u^{t} for u^{t+1} where A is the system matrix, solveLinearSystem(b,x)
    this->solveLinearSystem();
    //ierr = VecCopy(rightHandSide_, solution_); CHKERRV(ierr);  // for debugging just copy rhs (Vm_k^{*}) to solution (Vm_k^{i+1})

    // write the solution from the nested vector back to data
    int nSubVectors = 0;
    Vec *subVectors;
    ierr = VecNestGetSubVecs(solution_, &nSubVectors, &subVectors); CHKERRV(ierr);
    assert(nSubVectors == nCompartments_+1);

    for (int k = 0; k < nCompartments_; k++)
    {
      // copy the transmembrane potential from the subvector back to the subcellularStates vector, component 0 which is Vm
      dataMultidomain_.subcellularStates(k)->setValues(0, subVectors[k]);
    }

    VLOG(1) << "copy phi_e";

    // get phi_e
    dataMultidomain_.extraCellularPotential()->setValues(subVectors[nCompartments_]);

    LOG(DEBUG) << " Vm: ";
    dataMultidomain_.subcellularStates(0)->extractComponent(0, dataMultidomain_.transmembranePotential(0));
    LOG(DEBUG) << *dataMultidomain_.transmembranePotential(0);

    LOG(DEBUG) << " extraCellularPotential: ";
    LOG(DEBUG) << *dataMultidomain_.extraCellularPotential();

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

  LOG(INFO) << "Run potential flow simulation for fiber directions.";

  // solve potential flow Laplace problem
  finiteElementMethodPotentialFlow_.run();

  LOG(DEBUG) << "compute gradient field";

  // compute a gradient field from the solution of the potential flow
  dataMultidomain_.flowPotential()->setValues(*finiteElementMethodPotentialFlow_.data().solution());
  dataMultidomain_.flowPotential()->computeGradientField(dataMultidomain_.fibreDirection());

  LOG(DEBUG) << "flow potential: " << *dataMultidomain_.flowPotential();
  LOG(DEBUG) << "fiber direction: " << *dataMultidomain_.fibreDirection();

  // initialize the finite element class, from which only the stiffness matrix is needed
  finiteElementMethodDiffusion_.initialize(dataMultidomain_.fibreDirection());
  finiteElementMethodDiffusion_.initializeForImplicitTimeStepping(); // this performs extra initialization for implicit timestepping methods, i.e. it sets the inverse lumped mass matrix
  finiteElementMethodDiffusionTotal_.initialize(dataMultidomain_.fibreDirection(), true);

  initializeCellMLAdapters();

  // parse parameters
  PythonUtility::getOptionVector(this->specificSettings_, "am", nCompartments_, am_);
  PythonUtility::getOptionVector(this->specificSettings_, "cm", nCompartments_, cm_);
  LOG(DEBUG) << "Am: " << am_ << ", Cm: " << cm_;

  initializeCompartmentRelativeFactors();

  // initialize system matrix
  setSystemMatrix(this->timeStepWidth_);

  LOG(DEBUG) << "initialize linear solver";

  // initialize linear solver
  if (linearSolver_ == nullptr)
  {
    // retrieve linear solver
    linearSolver_ = this->context_.solverManager()->template solver<Solver::Linear>(
      this->specificSettings_, this->rankSubset_->mpiCommunicator());
  }

  LOG(DEBUG) << "set system matrix to linear solver";

  // set matrix used for linear solver and preconditioner to ksp context
  assert(this->linearSolver_->ksp());
  PetscErrorCode ierr;
  ierr = KSPSetOperators(*this->linearSolver_->ksp(), systemMatrix_, systemMatrix_); CHKERRV(ierr);

  // initialize rhs and solution vector
  subvectorsRightHandSide_.resize(nCompartments_+1);
  subvectorsSolution_.resize(nCompartments_+1);

  // set values for Vm in compartments
  for (int k = 0; k < nCompartments_; k++)
  {
    subvectorsRightHandSide_[k] = dataMultidomain_.ionicCurrent(k)->valuesGlobal();
    subvectorsSolution_[k] = dataMultidomain_.subcellularStates(k)->valuesGlobal(0);
  }

  // set values for phi_e
  subvectorsRightHandSide_[nCompartments_] = dataMultidomain_.zero()->valuesGlobal();
  subvectorsSolution_[nCompartments_] = dataMultidomain_.extraCellularPotential()->valuesGlobal();
  ierr = VecZeroEntries(subvectorsSolution_[nCompartments_]); CHKERRV(ierr);

  // create the nested vectors
  LOG(DEBUG) << "create nested vector";
  ierr = VecCreateNest(MPI_COMM_WORLD, nCompartments_+1, NULL, subvectorsRightHandSide_.data(), &rightHandSide_); CHKERRV(ierr);
  ierr = VecCreateNest(MPI_COMM_WORLD, nCompartments_+1, NULL, subvectorsSolution_.data(), &solution_); CHKERRV(ierr);

  LOG(DEBUG) << "initialization done";
  this->initialized_ = true;
}

template<typename FiniteElementMethodPotentialFlow,typename CellMLAdapterType,typename FiniteElementMethodDiffusion>
void MultidomainSolver<FiniteElementMethodPotentialFlow,CellMLAdapterType,FiniteElementMethodDiffusion>::
initializeCompartmentRelativeFactors()
{
  if (PythonUtility::hasKey(this->specificSettings_, "compartmentRelativeFactors"))
  {
    LOG(FATAL) << "compartmentRelativeFactors not yet implemented";
  }
  else
  {
    double fr = 1./nCompartments_;

    //initialize vectors with default values
    for (int k = 0; k < nCompartments_; k++)
    {
      dataMultidomain_.compartmentRelativeFactor(k)->setValues(fr);
    }
  }
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
    cellMLAdapters_[k].setInitialValues(dataMultidomain_.subcellularStates(k));
  }
}

template<typename FiniteElementMethodPotentialFlow,typename CellMLAdapterType,typename FiniteElementMethodDiffusion>
void MultidomainSolver<FiniteElementMethodPotentialFlow,CellMLAdapterType,FiniteElementMethodDiffusion>::
setSystemMatrix(double timeStepWidth)
{
  std::vector<Mat> submatrices(MathUtility::sqr(nCompartments_+1),NULL);

  LOG(TRACE) << "setSystemMatrix";

  // fill submatrices, empty submatrices may be NULL

  PetscErrorCode ierr;
  Mat stiffnessMatrix = finiteElementMethodDiffusion_.data().stiffnessMatrix()->valuesGlobal();
  Mat inverseLumpedMassMatrix = finiteElementMethodDiffusion_.data().inverseLumpedMassMatrix()->valuesGlobal();

  // set all submatrices
  for (int k = 0; k < nCompartments_; k++)
  {
    // right column matrix
    double prefactor = -this->timeStepWidth_ / (am_[k]*cm_[k]);

    VLOG(1) << "k=" << k << ", am: " << am_[k] << ", cm: " << cm_[k] << ", prefactor: " << prefactor;

    // create matrix as M^{-1}*K
    Mat matrixOnRightColumn;
    ierr = MatMatMult(inverseLumpedMassMatrix, stiffnessMatrix, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &matrixOnRightColumn);

    // scale matrix on right column with prefactor
    ierr = MatScale(matrixOnRightColumn, prefactor); CHKERRV(ierr);

    // copy right block matrix also to diagonal matrix
    Mat matrixOnDiagonalBlock;
    ierr = MatConvert(matrixOnRightColumn, MATSAME, MAT_INITIAL_MATRIX, &matrixOnDiagonalBlock); CHKERRV(ierr);

    // for debugging zero all entries
#ifdef MONODOMAIN
/**/    //ierr = MatZeroEntries(matrixOnRightColumn); CHKERRV(ierr);
#endif

    // set on right column of the system matrix
    submatrices[k*(nCompartments_+1) + (nCompartments_+1) - 1] = matrixOnRightColumn;

    VLOG(1) << "matrixOnRightColumn: " << PetscUtility::getStringMatrix(matrixOnRightColumn);
    VLOG(1) << "set at index " << k*(nCompartments_+1) + (nCompartments_+1) - 1;

    // ---
    // diagonal matrix
    // add identity
    ierr = MatShift(matrixOnDiagonalBlock, 1); CHKERRV(ierr);

    // set on diagonal
    submatrices[k*(nCompartments_+1) + k] = matrixOnDiagonalBlock;

    VLOG(1) << "matrixOnDiagonalBlock:" << PetscUtility::getStringMatrix(matrixOnDiagonalBlock);
    VLOG(1) << "set at index " << k*(nCompartments_+1) + k;

    // ---
    // bottom row matrices
    // create matrix as copy of stiffnessMatrix
    Mat matrixOnBottomRow;
    ierr = MatConvert(stiffnessMatrix, MATSAME, MAT_INITIAL_MATRIX, &matrixOnBottomRow); CHKERRV(ierr);

    // get the relative factor of the compartment
    Vec compartmentRelativeFactor = dataMultidomain_.compartmentRelativeFactor(k)->valuesGlobal();

    // MatDiagonalScale(A,l,NULL) computes A = diag(l)*A, this scales the rows of A with the values in l (each row with one entry of l)
    ierr = MatDiagonalScale(matrixOnBottomRow, compartmentRelativeFactor, NULL); CHKERRV(ierr);

    VLOG(1) << "matrixOnBottomRow: " << PetscUtility::getStringMatrix(matrixOnBottomRow);

    // for debugging zero all entries
#ifdef MONODOMAIN
/**/    ierr = MatZeroEntries(matrixOnBottomRow); CHKERRV(ierr);
#endif

    // set on bottom row of the system matrix
    submatrices[((nCompartments_+1) - 1)*(nCompartments_+1) + k] = matrixOnBottomRow;

    VLOG(1) << "set at index " << ((nCompartments_+1) - 1)*(nCompartments_+1) + k;
  }

  // set bottom right matrix
  Mat stiffnessMatrixBottomRight = finiteElementMethodDiffusionTotal_.data().stiffnessMatrix()->valuesGlobal();

  VLOG(1) << "stiffnessMatrixBottomRight:" << PetscUtility::getStringMatrix(stiffnessMatrixBottomRight);

  // for debugging set to identity
#ifdef MONODOMAIN
/**/  ierr = MatZeroEntries(stiffnessMatrixBottomRight); CHKERRV(ierr);
/**/  ierr = MatShift(stiffnessMatrixBottomRight,1); CHKERRV(ierr);
#endif

  // set on bottom right
  submatrices[(nCompartments_+1)*(nCompartments_+1)-1] = stiffnessMatrixBottomRight;
  VLOG(1) << "set at index " << (nCompartments_+1)*(nCompartments_+1)-1;

  // create nested matrix
  ierr = MatCreateNest(this->dataMultidomain_.functionSpace()->meshPartition()->mpiCommunicator(),
                       nCompartments_+1, NULL, nCompartments_+1, NULL, submatrices.data(), &this->systemMatrix_); CHKERRV(ierr);
}

template<typename FiniteElementMethodPotentialFlow,typename CellMLAdapterType,typename FiniteElementMethodDiffusion>
void MultidomainSolver<FiniteElementMethodPotentialFlow,CellMLAdapterType,FiniteElementMethodDiffusion>::
solveLinearSystem()
{
  PetscErrorCode ierr;


  VLOG(1) << "in solveLinearSystem";

  // configure that the initial value for the iterative solver is the value in solution, not zero
  ierr = KSPSetInitialGuessNonzero(*this->linearSolver_->ksp(), PETSC_TRUE); CHKERRV(ierr);

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
