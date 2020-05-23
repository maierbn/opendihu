#include "specialized_solver/multidomain_solver/multidomain_solver.h"

#include <Python.h>  // has to be the first included header
#include <iomanip>

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "data_management/specialized_solver/multidomain.h"
#include "specialized_solver/multidomain_solver/nested_mat_vec_utility.h"

//#define MONODOMAIN

namespace TimeSteppingScheme
{

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusion>
MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::
MultidomainSolver(DihuContext context) :
  TimeSteppingScheme(context["MultidomainSolver"]),
  dataMultidomain_(this->context_), finiteElementMethodPotentialFlow_(this->context_["PotentialFlow"]),
  finiteElementMethodDiffusion_(this->context_["Activation"]), finiteElementMethodDiffusionTotal_(this->context_["Activation"]),
  rankSubset_(DihuContext::partitionManager()->nextRankSubset()), nColumnSubmatricesSystemMatrix_(0)
{
  // get python config
  this->specificSettings_ = this->context_.getPythonConfig();

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_, this->specificSettings_);

  // parse number of motor units
  nCompartments_ = this->specificSettings_.getOptionInt("nCompartments", 1, PythonUtility::NonNegative);
  initialGuessNonzero_ = this->specificSettings_.getOptionBool("initialGuessNonzero", true);
  showLinearSolverOutput_ = this->specificSettings_.getOptionBool("showLinearSolverOutput", true);
  updateSystemMatrixEveryTimestep_ = this->specificSettings_.getOptionBool("updateSystemMatrixEveryTimestep", false);

  if (this->specificSettings_.hasKey("constructPreconditionerMatrix"))
  {
    LOG(ERROR) << this->specificSettings_ << " option \"constructPreconditionerMatrix\" has been renamed to \"useSymmetricPreconditionerMatrix\".";
  }
  useSymmetricPreconditionerMatrix_ = this->specificSettings_.getOptionBool("useSymmetricPreconditionerMatrix", true);

  // create finiteElement objects for diffusion in compartments
  finiteElementMethodDiffusionCompartment_.reserve(nCompartments_);
  for (int k = 0; k < nCompartments_; k++)
  {
    finiteElementMethodDiffusionCompartment_.emplace_back(this->context_["Activation"]);
  }

  singleSystemMatrix_ = PETSC_NULL;
  singleSolution_ = PETSC_NULL;
  singleRightHandSide_ = PETSC_NULL;
  singlePreconditionerMatrix_ = PETSC_NULL;
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusion>
void MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::
advanceTimeSpan()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  // compute timestep width
  double timeSpan = this->endTime_ - this->startTime_;

  LOG(DEBUG) << "MultidomainSolver::advanceTimeSpan, timeSpan=" << timeSpan<< ", timeStepWidth=" << this->timeStepWidth_
    << " n steps: " << this->numberTimeSteps_;

  // loop over time steps
  double currentTime = this->startTime_;

  // loop over time steps
  for (int timeStepNo = 0; timeStepNo < this->numberTimeSteps_;)
  {
    if (timeStepNo % this->timeStepOutputInterval_ == 0 && (this->timeStepOutputInterval_ <= 10 || timeStepNo > 0))  // show first timestep only if timeStepOutputInterval is <= 10
    {
      LOG(INFO) << "Multidomain diffusion, timestep " << timeStepNo << "/" << this->numberTimeSteps_<< ", t=" << currentTime
        << " (linear solver iterations: " << lastNumberOfIterations_ << ")";
    }

    LOG(DEBUG) << " Vm: ";
    //dataMultidomain_.subcellularStates(0)->extractComponent(0, dataMultidomain_.transmembranePotential(0));
    LOG(DEBUG) << *dataMultidomain_.transmembranePotential(0);

    if (fabs(this->timeStepWidthOfSystemMatrix_ - this->timeStepWidth_) / this->timeStepWidth_ > 1e-4)
    {
      //LOG(WARNING) << "In multidomain solver, timestep width changed from " << this->timeStepWidthOfSystemMatrix_ << " to " << timeStepWidth_
      //  << " (relative: " << std::showpos << 100*(this->timeStepWidthOfSystemMatrix_ - this->timeStepWidth_) / this->timeStepWidth_ << std::noshowpos << "%), need to recreate system matrix.";

      this->timeStepWidthOfSystemMatrix_ = this->timeStepWidth_;
      setSystemMatrixSubmatrices(this->timeStepWidthOfSystemMatrix_);
      createSystemMatrixFromSubmatrices();
    }
    else if (this->updateSystemMatrixEveryTimestep_ && timeStepNo == 0)
    {
      updateSystemMatrix();
    }

    // advance simulation time
    timeStepNo++;
    currentTime = this->startTime_ + double(timeStepNo) / this->numberTimeSteps_ * timeSpan;

    // advance diffusion
    VLOG(1) << "---- diffusion term";

    // solve A*u^{t+1} = u^{t} for u^{t+1} where A is the system matrix, solveLinearSystem(b,x)
    this->solveLinearSystem();

    LOG(DEBUG) << " Vm[k=0]: ";
    //dataMultidomain_.subcellularStates(0)->extractComponent(0, dataMultidomain_.transmembranePotential(0));
    LOG(DEBUG) << *dataMultidomain_.transmembranePotentialSolution(0);

    LOG(DEBUG) << " extraCellularPotential: " << PetscUtility::getStringVector(subvectorsSolution_[nCompartments_]);
    LOG(DEBUG) << *dataMultidomain_.extraCellularPotential();

    // stop duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::stop(this->durationLogKey_);

    // write current output values
    callOutputWriter(timeStepNo, currentTime);

    // start duration measurement
    if (this->durationLogKey_ != "")
      Control::PerformanceMeasurement::start(this->durationLogKey_);
  }

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);

}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusion>
void MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::
run()
{
  // initialize everything
  initialize();

  this->advanceTimeSpan();
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusion>
void MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::
callOutputWriter(int timeStepNo, double currentTime)
{
  // write current output values
  this->outputWriterManager_.writeOutput(this->dataMultidomain_, timeStepNo, currentTime);
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusion>
void MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::
initialize()
{
  if (this->initialized_)
    return;

  initializeObjects();              // this is also called by MultidomainWithFatSolver
  initializeMatricesAndVectors();   // this is not called by MultidomainWithFatSolver

  // write initial meshes
  callOutputWriter(0, 0.0);

  LOG(DEBUG) << "initialization done";
  this->initialized_ = true;
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusion>
void MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::
initializeObjects()
{
  //! initialize everything except the matrices and vectors, this will also be called by the inherited class
  // call initialize of the timestepping scheme to set timestep width
  TimeSteppingScheme::initialize();

  LOG(DEBUG) << "initialize multidomain_solver, " << nCompartments_ << " compartments";
  assert(this->specificSettings_.pyObject());

  // add this solver to the solvers diagram
  DihuContext::solverStructureVisualizer()->addSolver("MultidomainSolver");

  // indicate in solverStructureVisualizer that now a child solver will be initialized
  DihuContext::solverStructureVisualizer()->beginChild("PotentialFlow");

  // initialize the potential flow finite element method, this also creates the function space
  finiteElementMethodPotentialFlow_.initialize();

  // indicate in solverStructureVisualizer that the child solver initialization is done
  DihuContext::solverStructureVisualizer()->endChild();

  // initialize the data object
  dataMultidomain_.setFunctionSpace(finiteElementMethodPotentialFlow_.functionSpace());
  dataMultidomain_.initialize(nCompartments_);

  DihuContext::solverStructureVisualizer()->setOutputConnectorData(getOutputConnectorData());

  LOG(INFO) << "Run potential flow simulation for fiber directions.";

  // solve potential flow Laplace problem
  finiteElementMethodPotentialFlow_.run();

  LOG(DEBUG) << "compute gradient field";

  // compute a gradient field from the solution of the potential flow
  dataMultidomain_.flowPotential()->setValues(*finiteElementMethodPotentialFlow_.data().solution());
  dataMultidomain_.flowPotential()->computeGradientField(dataMultidomain_.fiberDirection());

  VLOG(1) << "flow potential: " << *dataMultidomain_.flowPotential();
  VLOG(1) << "fiber direction: " << *dataMultidomain_.fiberDirection();

  initializeCompartmentRelativeFactors();

  // indicate in solverStructureVisualizer that now a child solver will be initialized
  DihuContext::solverStructureVisualizer()->beginChild("Activation");

  // initialize the finite element class, from which only the stiffness matrix is needed
  // diffusion object without prefactor, for normal diffusion (2nd multidomain eq.)
  finiteElementMethodDiffusion_.initialize(dataMultidomain_.fiberDirection(), nullptr);
  finiteElementMethodDiffusion_.initializeForImplicitTimeStepping(); // this performs extra initialization for implicit timestepping methods, i.e. it sets the inverse lumped mass matrix

  // indicate in solverStructureVisualizer that the child solver initialization is done
  DihuContext::solverStructureVisualizer()->endChild();
  // do not log the compartment finite element method objects
  DihuContext::solverStructureVisualizer()->disable();

  // system to be solved:
  //
  // [A^1_Vm,Vm   |            |             | B^1_Vm,phie]   [ V^1_m^(i+1) ]   [V^1_m^(i)]
  // [            | A^2_Vm,Vm  |             | B^2_Vm,phie] * [ V^2_m^(i+1) ] = [V^2_m^(i)]
  // [   ...      |            |  A^M_Vm,Vm  | B^M_Vm,phie]   [ V^M_m^(i+1) ]   [V^M_m^(i)]
  // [B^1_phie,Vm |B^2_phie,Vm | B^M_phie,Vm | B_phie,phie]   [ phi_e^(i+1) ]   [0        ]

  // diffusion objects with spatially varying prefactors (f_r), needed for the bottom row of the matrix eq. or the 1st multidomain eq.
  for (int k = 0; k < nCompartments_; k++)
  {
    finiteElementMethodDiffusionCompartment_[k].initialize(dataMultidomain_.fiberDirection(), dataMultidomain_.compartmentRelativeFactor(k));
    finiteElementMethodDiffusionCompartment_[k].initializeForImplicitTimeStepping(); // this performs extra initialization for implicit timestepping methods, i.e. it sets the inverse lumped mass matrix
  }

  finiteElementMethodDiffusionTotal_.initialize(dataMultidomain_.fiberDirection(), dataMultidomain_.relativeFactorTotal(), true);

  DihuContext::solverStructureVisualizer()->enable();

  // parse parameters
  this->specificSettings_.getOptionVector("am", nCompartments_, am_);
  this->specificSettings_.getOptionVector("cm", nCompartments_, cm_);
  LOG(DEBUG) << "Am: " << am_ << ", Cm: " << cm_;

  // initialize linear solver
  LOG(DEBUG) << "initialize linear solver";

  if (this->linearSolver_ == nullptr)
  {
    // create or get linear solver object
    this->linearSolver_ = this->context_.solverManager()->template solver<Solver::Linear>(
      this->specificSettings_, this->rankSubset_->mpiCommunicator());

    // initialize the alternative linear solver that is used when thet linearSolver_ diverges
    if (this->specificSettings_.hasKey("alternativeSolverName"))
    {
      this->alternativeLinearSolver_ = this->context_.solverManager()->template solver<Solver::Linear>(
        this->specificSettings_, this->rankSubset_->mpiCommunicator(), "alternativeSolverName");
    }
  }
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusion>
void MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::
initializeMatricesAndVectors()
{
  LOG(DEBUG) << "initialize linear solver";

  // initialize system matrix
  updateSystemMatrix(this->timeStepWidth_, false);

  LOG(DEBUG) << "set system matrix to linear solver";

  // set matrix used for linear solver and preconditioner to ksp context
  assert(this->linearSolver_->ksp());
  PetscErrorCode ierr;

  PC pc;
  ierr = KSPGetPC(*linearSolver_->ksp(), &pc); CHKERRV(ierr);

  // set block information for block jacobi preconditioner
  // check, if block jacobi preconditioner is selected
  PetscBool useBlockJacobiPreconditioner;
  PetscObjectTypeCompare((PetscObject)pc, PCBJACOBI, &useBlockJacobiPreconditioner);
  if (useBlockJacobiPreconditioner)
  {
    // PCBJacobiSetTotalBlocks(PC pc, PetscInt nBlocks, const PetscInt lengthsOfBlocks[])
    PetscInt nBlocks = nColumnSubmatricesSystemMatrix_;

    // set sizes of all blocks to the number of dofs in the muscle domain
    std::vector<PetscInt> lengthsOfBlocks(nBlocks, dataMultidomain_.functionSpace()->nDofsGlobal());
    ierr = PCBJacobiSetTotalBlocks(pc, nColumnSubmatricesSystemMatrix_, lengthsOfBlocks.data()); CHKERRV(ierr);
  }

  // set the local node positions for the preconditioner
  int nDofsPerNode = dataMultidomain_.functionSpace()->nDofsPerNode();
  int nNodesLocal = dataMultidomain_.functionSpace()->nNodesLocalWithoutGhosts();

  std::vector<double> nodePositionCoordinatesForPreconditioner;
  nodePositionCoordinatesForPreconditioner.reserve(3*nNodesLocal);

  // loop over muscle nodes and add their node positions
  for (dof_no_t dofNoLocal = 0; dofNoLocal < nNodesLocal*nDofsPerNode; dofNoLocal++)
  {
    Vec3 nodePosition = dataMultidomain_.functionSpace()->getGeometry(dofNoLocal);

    // add the coordinates
    for (int i = 0; i < 3; i++)
      nodePositionCoordinatesForPreconditioner.push_back(nodePosition[i]);
  }

  LOG(DEBUG) << "set coordinates to preconditioner, " << nodePositionCoordinatesForPreconditioner.size() << " node coordinates";
  ierr = PCSetCoordinates(pc, 3, nodePositionCoordinatesForPreconditioner.size(), nodePositionCoordinatesForPreconditioner.data()); CHKERRV(ierr);

  // set the nullspace of the matrix
  // as we have Neumann boundary conditions, constant functions are in the nullspace of the matrix
  MatNullSpace nullSpace;
  ierr = MatNullSpaceCreate(data().functionSpace()->meshPartition()->mpiCommunicator(), PETSC_TRUE, 0, PETSC_NULL, &nullSpace); CHKERRV(ierr);
  ierr = MatSetNullSpace(singleSystemMatrix_, nullSpace); CHKERRV(ierr);
  ierr = MatSetNearNullSpace(singleSystemMatrix_, nullSpace); CHKERRV(ierr); // for multigrid methods
  //ierr = MatNullSpaceDestroy(&nullSpace); CHKERRV(ierr);

  ierr = KSPSetOperators(*this->linearSolver_->ksp(), singleSystemMatrix_, singlePreconditionerMatrix_); CHKERRV(ierr);

  // initialize rhs and solution vector
  subvectorsRightHandSide_.resize(nCompartments_+1);
  subvectorsSolution_.resize(nCompartments_+1);

  // set values for Vm in compartments
  for (int k = 0; k < nCompartments_; k++)
  {
    subvectorsRightHandSide_[k] = dataMultidomain_.transmembranePotential(k)->valuesGlobal();     // this is for V_mk^(i)
    subvectorsSolution_[k] = dataMultidomain_.transmembranePotentialSolution(k)->valuesGlobal(0); // this is for V_mk^(i+1)
  }

  // set values for phi_e
  subvectorsRightHandSide_[nCompartments_] = dataMultidomain_.zero()->valuesGlobal();
  subvectorsSolution_[nCompartments_] = dataMultidomain_.extraCellularPotential()->valuesGlobal();
  ierr = VecZeroEntries(subvectorsSolution_[nCompartments_]); CHKERRV(ierr);

  // create the nested vectors
  LOG(DEBUG) << "create nested vector";
  ierr = VecCreateNest(this->rankSubset_->mpiCommunicator(), nCompartments_+1, NULL, subvectorsRightHandSide_.data(), &nestedRightHandSide_); CHKERRV(ierr);
  ierr = VecCreateNest(this->rankSubset_->mpiCommunicator(), nCompartments_+1, NULL, subvectorsSolution_.data(), &nestedSolution_); CHKERRV(ierr);

  // copy the values from a nested Petsc Vec to a single Vec that contains all entries
  NestedMatVecUtility::createVecFromNestedVec(nestedRightHandSide_, singleRightHandSide_, data().functionSpace()->meshPartition()->rankSubset());
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusion>
void MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::
initializeCompartmentRelativeFactors()
{
  // parse relative factors f_r for compartments
  bool inputIsGlobal = this->specificSettings_.getOptionBool("inputIsGlobal", true);

  std::vector<PythonConfig> compartmentFields;
  this->specificSettings_.getOptionVector("compartmentRelativeFactors", compartmentFields);
  if (compartmentFields.size() < nCompartments_)
  {
    LOG(FATAL) << "Only " << compartmentFields.size() << " relative factors specified under \"compartmentRelativeFactors\". "
      << "Number of compartments is " << nCompartments_ << ".";
  }


  for (int k = 0; k < nCompartments_; k++)
  {
    std::vector<double> values = PythonUtility::convertFromPython<std::vector<double>>::get(compartmentFields[k].pyObject());

    // if parsed node positions in vector localNodePositions_ actually contains global node positions, extract local positions
    if (inputIsGlobal)
    {
      if (values.size() != dataMultidomain_.functionSpace()->nDofsGlobal())
      {
        LOG(FATAL) << "In MultidomainSolver, \"compartmentRelativeFactors\" for compartment " << k << " of " << nCompartments_ << " contains "
          << values.size() << " entries, but the mesh \"" << dataMultidomain_.functionSpace()->meshName() << "\" has "
          << dataMultidomain_.functionSpace()->nDofsGlobal() << " global dofs and \"inputIsGlobal\" is True.\n"
          << "Depending on how you pass the values in the python settings, maybe delete and recreate the \"compartments_relative_factors*\" files?\n"
          << "(Note this is a guess, the C++ code does not know which example you're in or what you're doing in the python scripts.)";
      }

      dataMultidomain_.functionSpace()->meshPartition()->extractLocalDofsWithoutGhosts(values);
    }

    if (values.size() < dataMultidomain_.compartmentRelativeFactor(k)->nDofsLocalWithoutGhosts())
    {
      LOG(FATAL) << "In MultidomainSolver, \"compartmentRelativeFactors\" for compartment " << k << " of " << nCompartments_
        << " contains only " << values.size() << " entries, but the mesh \"" << dataMultidomain_.compartmentRelativeFactor(k)->functionSpace()->meshName() << "\""
        << " has " << dataMultidomain_.compartmentRelativeFactor(k)->nDofsLocalWithoutGhosts() << " local dofs.\n"
        << "Depending on how you pass the values in the python settings, maybe delete and recreate the \"compartments_relative_factors*\" files?\n"
        << "(Note this is a guess, the C++ code does not know which example you're in or what you're doing in the python scripts.)";
    }

    dataMultidomain_.compartmentRelativeFactor(k)->setValuesWithoutGhosts(values);
  }

  // compute relative Factor total as sum f_r
  dataMultidomain_.relativeFactorTotal()->zeroEntries();
  PetscErrorCode ierr;
  for (int k = 0; k < nCompartments_; k++)
  {
    ierr = VecAXPY(dataMultidomain_.relativeFactorTotal()->valuesGlobal(), 1.0, dataMultidomain_.compartmentRelativeFactor(k)->valuesGlobal()); CHKERRV(ierr);
  }


  for (int k = 0; k < nCompartments_; k++)
  {
    LOG(DEBUG) << "compartmentRelativeFactor(k=" << k << "): " << *dataMultidomain_.compartmentRelativeFactor(k);
  }
  LOG(DEBUG) << "relativeFactorTotal: " << *dataMultidomain_.relativeFactorTotal();

}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusion>
void MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::
setSystemMatrixSubmatrices(double timeStepWidth)
{
  //LOG(INFO) << "dt of system matrix: " << timeStepWidth;
  // initialize the number of rows, this will already be set when the method is called for the inherited class MultidomainWithFatSolver
  if (nColumnSubmatricesSystemMatrix_ == 0)
    nColumnSubmatricesSystemMatrix_ = nCompartments_+1;

  this->submatricesSystemMatrix_.resize(MathUtility::sqr(nColumnSubmatricesSystemMatrix_),NULL);

  LOG(TRACE) << "setSystemMatrix";

  // The system to be solved here is
  //
  // [ -dt/(a_mk*c_mk)*M^{-1}*K_ik + I   ...   -dt/(a_mk*c_mk)*M^{-1}*K   ] [V_mk^(i+1)  ]   [V_mk^(i)]
  // [  ...                                     ...                       ]*[...         ] = [       ]
  // [ f_rk * K_ik                       ...    K_ei                      ] [phi_e^(i+1) ]   [       ]
  //
  // V_mk^(i) is computed by the 0D part.
  // The two output connection slots are V_mk^(i) and V_mk^(i+1),
  // where V_mk^(i) is the input and should be connected to the output of the reaction term.
  // V_mk^(i+1) is the output and should be connected to the input of the reaction term.

  // fill submatricesSystemMatrix_, empty submatricesSystemMatrix_ may be NULL
  // stiffnessMatrix and inverse lumped mass matrix without prefactor
  Mat stiffnessMatrix = finiteElementMethodDiffusion_.data().stiffnessMatrix()->valuesGlobal();
  Mat inverseLumpedMassMatrix = finiteElementMethodDiffusion_.data().inverseLumpedMassMatrix()->valuesGlobal();

  PetscErrorCode ierr;
  // set all submatricesSystemMatrix_
  for (int k = 0; k < nCompartments_; k++)
  {

    // right column matrix
    double prefactor = -timeStepWidth / (am_[k]*cm_[k]);

    VLOG(2) << "k=" << k << ", am: " << am_[k] << ", cm: " << cm_[k] << ", prefactor: " << prefactor;

    // create matrix as M^{-1}*K
    Mat matrixOnRightColumn;
    ierr = MatMatMult(inverseLumpedMassMatrix, stiffnessMatrix, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &matrixOnRightColumn); CHKERRV(ierr);

    // scale matrix on right column with prefactor
    ierr = MatScale(matrixOnRightColumn, prefactor); CHKERRV(ierr);

    // copy right block matrix also to diagonal matrix
    Mat matrixOnDiagonalBlock;
    ierr = MatConvert(matrixOnRightColumn, MATSAME, MAT_INITIAL_MATRIX, &matrixOnDiagonalBlock); CHKERRV(ierr);

    // for debugging zero all entries
#ifdef MONODOMAIN
/**/    ierr = MatZeroEntries(matrixOnRightColumn); CHKERRV(ierr);
#endif

    // set on right column of the system matrix
    submatricesSystemMatrix_[k*nColumnSubmatricesSystemMatrix_ + (nCompartments_+1) - 1] = matrixOnRightColumn;

    if (VLOG_IS_ON(2))
    {
      VLOG(2) << "matrixOnRightColumn: " << PetscUtility::getStringMatrix(matrixOnRightColumn);
      VLOG(2) << "set at index " << k*nColumnSubmatricesSystemMatrix_ + (nCompartments_+1) - 1;
    }

    // ---
    // diagonal matrix
    // add identity
    ierr = MatShift(matrixOnDiagonalBlock, 1); CHKERRV(ierr);

    // set on diagonal
    submatricesSystemMatrix_[k*nColumnSubmatricesSystemMatrix_ + k] = matrixOnDiagonalBlock;

    if (VLOG_IS_ON(2))
    {
      VLOG(2) << "matrixOnDiagonalBlock:" << PetscUtility::getStringMatrix(matrixOnDiagonalBlock);
      VLOG(2) << "set at index " << k*nColumnSubmatricesSystemMatrix_ + k;
    }

    // ---
    // bottom row matrices

    // stiffnessMatrixWithPrefactor is f_k*K
    Mat stiffnessMatrixWithPrefactor = finiteElementMethodDiffusionCompartment_[k].data().stiffnessMatrix()->valuesGlobal();

    // create matrix as copy of stiffnessMatrix
    Mat matrixOnBottomRow;
    ierr = MatConvert(stiffnessMatrixWithPrefactor, MATSAME, MAT_INITIAL_MATRIX, &matrixOnBottomRow); CHKERRV(ierr);

#if 0
    // debugging test, gives slightly different results due to approximation of test
    Mat test;
    ierr = MatConvert(stiffnessMatrix, MATSAME, MAT_INITIAL_MATRIX, &test); CHKERRV(ierr);

    // get the relative factor of the compartment
    Vec compartmentRelativeFactor = dataMultidomain_.compartmentRelativeFactor(k)->valuesGlobal();

    int i = 0;double v = 0;
    MatGetValues(test,1,&i,1,&i,&v);
    LOG(DEBUG) << "compartmentRelativeFactor: " << PetscUtility::getStringVector(compartmentRelativeFactor);
    LOG(DEBUG) << "test[0][0] = " << v;

    // MatDiagonalScale(A,l,NULL) computes A = diag(l)*A, this scales the rows of A with the values in l (each row with one entry of l)
    ierr = MatDiagonalScale(test, compartmentRelativeFactor, NULL); CHKERRV(ierr);

    MatGetValues(test,1,&i,1,&i,&v);
    LOG(DEBUG) << "MatdiagonalScale";
    LOG(DEBUG) << "test[0][0] = " << v;

    // test if "test" and matrixOnBottomRow yield the same matrices
    PetscInt nRows, nColumns;
    ierr = MatGetSize(matrixOnBottomRow,&nRows,&nColumns); CHKERRV(ierr);
    double max_diff = 0;
    for (int column = 0; column < nColumns; column++)
    {
      for (int row = 0; row < nRows; row++)
      {
        double value1, value2;
        ierr = MatGetValues(matrixOnBottomRow,1,&row,1,&column,&value1); CHKERRV(ierr);
        ierr = MatGetValues(test,1,&row,1,&column,&value2); CHKERRV(ierr);
        double difference = fabs(value1-value2);
        max_diff = std::max(max_diff, difference);
        if (difference > 1e-5)
        {
          LOG(DEBUG) << "bottom matrix entry row=" << row << ", column=" << column << " is different: " << value1 << "," << value2 << ", diff: " << difference;

        }
      }
    }
    LOG(FATAL) << "max_diff:" << max_diff;
#endif
    if (VLOG_IS_ON(2))
    {
      VLOG(2) << "matrixOnBottomRow: " << PetscUtility::getStringMatrix(matrixOnBottomRow);
    }

    // for debugging zero all entries
#ifdef MONODOMAIN
/**/    ierr = MatZeroEntries(matrixOnBottomRow); CHKERRV(ierr);
#endif

    // set on bottom row of the system matrix
    submatricesSystemMatrix_[((nCompartments_+1) - 1)*nColumnSubmatricesSystemMatrix_ + k] = matrixOnBottomRow;

    VLOG(2) << "set at index " << ((nCompartments_+1) - 1)*nColumnSubmatricesSystemMatrix_ + k;
  }

  // set bottom right matrix
  Mat stiffnessMatrixBottomRight = finiteElementMethodDiffusionTotal_.data().stiffnessMatrix()->valuesGlobal();

  if (VLOG_IS_ON(2))
  {
    VLOG(2) << "stiffnessMatrixBottomRight:" << PetscUtility::getStringMatrix(stiffnessMatrixBottomRight);
  }

  // for debugging set to identity
#ifdef MONODOMAIN
/**/  ierr = MatZeroEntries(stiffnessMatrixBottomRight); CHKERRV(ierr);
/**/  ierr = MatShift(stiffnessMatrixBottomRight,1); CHKERRV(ierr);
#endif

  // set on bottom right
  submatricesSystemMatrix_[((nCompartments_+1) - 1)*nColumnSubmatricesSystemMatrix_ + (nCompartments_+1)-1] = stiffnessMatrixBottomRight;
  VLOG(2) << "set at index " << ((nCompartments_+1) - 1)*nColumnSubmatricesSystemMatrix_ + (nCompartments_+1)-1;

}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusion>
void MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::
createSystemMatrixFromSubmatrices()
{
  PetscErrorCode ierr;
  assert(submatricesSystemMatrix_.size() == nColumnSubmatricesSystemMatrix_*nColumnSubmatricesSystemMatrix_);

#ifndef NDEBUG
  LOG(DEBUG) << "nested matrix with " << nColumnSubmatricesSystemMatrix_ << "x" << nColumnSubmatricesSystemMatrix_ << " submatrices, nCompartments_=" << nCompartments_;

  // output dimensions of submatrices for debugging
  for (int rowNo = 0; rowNo < nColumnSubmatricesSystemMatrix_; rowNo++)
  {
    for (int columnNo = 0; columnNo < nColumnSubmatricesSystemMatrix_; columnNo++)
    {
      Mat subMatrix = submatricesSystemMatrix_[rowNo*nColumnSubmatricesSystemMatrix_ + columnNo];

      if (!subMatrix)
      {
        LOG(DEBUG) << "submatrix (" << rowNo << "," << columnNo << ") is empty (NULL)";
      }
      else
      {
        PetscInt nRows, nColumns;
        ierr = MatGetSize(subMatrix, &nRows, &nColumns); CHKERRV(ierr);
        std::string name;
        char *cName;
        ierr = PetscObjectGetName((PetscObject)subMatrix, (const char **)&cName); CHKERRV(ierr);
        name = cName;

        LOG(DEBUG) << "submatrix (" << rowNo << "," << columnNo << ") is \"" << name << "\" (" << nRows << "x" << nColumns << ")";
      }
    }
  }
#endif

  // create nested matrix
  ierr = MatCreateNest(this->rankSubset_->mpiCommunicator(),
                       nColumnSubmatricesSystemMatrix_, NULL, nColumnSubmatricesSystemMatrix_, NULL, submatricesSystemMatrix_.data(), &nestedSystemMatrix_); CHKERRV(ierr);

  // create a single Mat object from the nested Mat
  NestedMatVecUtility::createMatFromNestedMat(nestedSystemMatrix_, singleSystemMatrix_, data().functionSpace()->meshPartition()->rankSubset());

  if (useSymmetricPreconditionerMatrix_)
  {
    this->submatricesPreconditionerMatrix_ = submatricesSystemMatrix_;

    // set offdiagonal matrices to NULL
    for (int rowNo = 0; rowNo < nColumnSubmatricesSystemMatrix_; rowNo++)
    {
      for (int columnNo = 0; columnNo < nColumnSubmatricesSystemMatrix_; columnNo++)
      {
        int index = rowNo*nColumnSubmatricesSystemMatrix_ + columnNo;

        bool isOnDiagonal = columnNo == rowNo;
        bool isSymmetricOffDiagonal = (rowNo == nCompartments_+1 && columnNo == nCompartments_) || (rowNo == nCompartments_ && columnNo == nCompartments_+1);

        if (!isOnDiagonal && !isSymmetricOffDiagonal)
        {
          submatricesPreconditionerMatrix_[index] = NULL;
        }
      }
    }
#ifndef NDEBUG
  LOG(DEBUG) << "preconditioner: nested matrix with " << nColumnSubmatricesSystemMatrix_ << "x" << nColumnSubmatricesSystemMatrix_ << " submatrices, nCompartments_=" << nCompartments_;

  // output dimensions of submatrices for debugging
  for (int rowNo = 0; rowNo < nColumnSubmatricesSystemMatrix_; rowNo++)
  {
    for (int columnNo = 0; columnNo < nColumnSubmatricesSystemMatrix_; columnNo++)
    {
      Mat subMatrix = submatricesPreconditionerMatrix_[rowNo*nColumnSubmatricesSystemMatrix_ + columnNo];

      if (!subMatrix)
      {
        LOG(DEBUG) << "preconditioner submatrix (" << rowNo << "," << columnNo << ") is empty (NULL)";
      }
      else
      {
        PetscInt nRows, nColumns;
        ierr = MatGetSize(subMatrix, &nRows, &nColumns); CHKERRV(ierr);
        std::string name;
        char *cName;
        ierr = PetscObjectGetName((PetscObject)subMatrix, (const char **)&cName); CHKERRV(ierr);
        name = cName;

        LOG(DEBUG) << "preconditioner submatrix (" << rowNo << "," << columnNo << ") is \"" << name << "\" (" << nRows << "x" << nColumns << ")";
      }
    }
  }
#endif


    Mat nestedPreconditionerMatrix;

    // create nested matrix
    ierr = MatCreateNest(this->rankSubset_->mpiCommunicator(),
                        nColumnSubmatricesSystemMatrix_, NULL, nColumnSubmatricesSystemMatrix_, NULL, submatricesPreconditionerMatrix_.data(), &nestedPreconditionerMatrix); CHKERRV(ierr);

    // create a single Mat object from the nested Mat
    NestedMatVecUtility::createMatFromNestedMat(nestedPreconditionerMatrix, singlePreconditionerMatrix_, data().functionSpace()->meshPartition()->rankSubset());
  }
  else
  {
    singlePreconditionerMatrix_ = singleSystemMatrix_;
  }
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusion>
void MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::
updateSystemMatrix()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_+std::string("_reassemble"));

  // assemble stiffness and mass matrices again
  this->finiteElementMethodDiffusion_.setStiffnessMatrix();
  this->finiteElementMethodDiffusion_.setMassMatrix();
  this->finiteElementMethodDiffusion_.setInverseLumpedMassMatrix();

  this->finiteElementMethodDiffusionTotal_.setStiffnessMatrix();

  // compute new entries for submatrices, except B,C,D and E
  setSystemMatrixSubmatrices(this->timeStepWidthOfSystemMatrix_);

  // create the system matrix again
  createSystemMatrixFromSubmatrices();

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_+std::string("_reassemble"));
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusion>
void MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::
updateSystemMatrix(double timeStepWidth, bool enableWarning)
{
  if (fabs(this->timeStepWidthOfSystemMatrix_ - timeStepWidth) / timeStepWidth > 1e-4)
  {
    if (enableWarning)
    {
      //LOG(WARNING) << "In multidomain solver, timestep width changed from " << this->timeStepWidthOfSystemMatrix_ << " to " << timeStepWidth_
      //  << " (relative: " << std::showpos << 100*(this->timeStepWidthOfSystemMatrix_ - timeStepWidth) / timeStepWidth << std::noshowpos << "%), need to recreate system matrix.";
    }

    // store time step width of current system matrix
    this->timeStepWidthOfSystemMatrix_ = timeStepWidth;

    // initialize the sub matrices of the nested system matrix
    setSystemMatrixSubmatrices(this->timeStepWidthOfSystemMatrix_);

    // create the nested and the single system matrix
    createSystemMatrixFromSubmatrices();
  }
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusion>
void MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::
solveLinearSystem()
{
  VLOG(1) << "in solveLinearSystem";

  // configure that the initial value for the iterative solver is the value in solution, not zero
  PetscErrorCode ierr;
  if (initialGuessNonzero_)
  {
    LOG(DEBUG) << "set initial guess nonzero";
    ierr = KSPSetInitialGuessNonzero(*this->linearSolver_->ksp(), PETSC_TRUE); CHKERRV(ierr);
  }

  // copy the values from a nested Petsc Vec to a single Vec that contains all entries
  NestedMatVecUtility::createVecFromNestedVec(nestedRightHandSide_, singleRightHandSide_, data().functionSpace()->meshPartition()->rankSubset());

  // solve the linear system
  // this can be done using the nested Vecs and nested Mat (nestedSolution_, nestedRightHandSide_, nestedSystemMatrix_),
  // or the single Vecs and Mats that contain all values directly  (singleSolution_, singleRightHandSide_, singleSystemMatrix_)

  bool hasSolverConverged = false;

  // try up to three times to solve the system
  for (int solveNo = 0; solveNo < 5; solveNo++)
  {
    // copy the values from a nested Petsc Vec to a single Vec that contains all entries
    NestedMatVecUtility::createVecFromNestedVec(nestedSolution_, singleSolution_, data().functionSpace()->meshPartition()->rankSubset());

    if (showLinearSolverOutput_)
    {
      hasSolverConverged = this->linearSolver_->solve(singleRightHandSide_, singleSolution_, "Linear system of multidomain problem solved");
    }
    else
    {
      // solve without showing output
      hasSolverConverged = this->linearSolver_->solve(singleRightHandSide_, singleSolution_);
    }
    if (hasSolverConverged)
    {
      break;
    }
    else
    {
      LOG(WARNING) << "Solver has not converged, try again " << solveNo << "/3";
    }
  }

  // store the last number of iterations
  lastNumberOfIterations_ = this->linearSolver_->lastNumberOfIterations();

  // copy the values back from a single Vec that contains all entries to a nested Petsc Vec
  NestedMatVecUtility::fillNestedVec(singleSolution_, nestedSolution_);
}

template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusion>
typename MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::Data &MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::
data()
{
  return dataMultidomain_;
}

//! get the data that will be transferred in the operator splitting to the other term of the splitting
//! the transfer is done by the output_connector_data_transfer class
template<typename FiniteElementMethodPotentialFlow,typename FiniteElementMethodDiffusion>
std::shared_ptr<typename MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::OutputConnectorDataType>
MultidomainSolver<FiniteElementMethodPotentialFlow,FiniteElementMethodDiffusion>::
getOutputConnectorData()
{
  LOG(DEBUG) << "getOutputConnectorData, size of Vm vector: " << this->dataMultidomain_.transmembranePotential().size();

  return dataMultidomain_.getOutputConnectorData();
}

} // namespace TimeSteppingScheme
