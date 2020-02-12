#include "specialized_solver/solid_mechanics/hyperelasticity/hyperelasticity_solver.h"

#include <Python.h>  // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "solver/solver_manager.h"
#include "data_management/specialized_solver/multidomain.h"
#include "control/diagnostic_tool/performance_measurement.h"
#include "control/diagnostic_tool/solver_structure_visualizer.h"

namespace SpatialDiscretization
{

template<typename Term,int nDisplacementComponents>
HyperelasticitySolver<Term,nDisplacementComponents>::
HyperelasticitySolver(DihuContext context, std::string settingsKey) :
  context_(context[settingsKey]), data_(context_), pressureDataCopy_(context_), initialized_(false),
  endTime_(0), lastNorm_(0), secondLastNorm_(0)
{
  // get python config
  this->specificSettings_ = this->context_.getPythonConfig();

  // read in the durationLogKey
  if (specificSettings_.hasKey("durationLogKey"))
  {
    this->durationLogKey_ = specificSettings_.getOptionString("durationLogKey", "");
  }

  // parse options concerning jacobian
  useAnalyticJacobian_ = this->specificSettings_.getOptionBool("useAnalyticJacobian", true);
  useNumericJacobian_ = this->specificSettings_.getOptionBool("useNumericJacobian", true);

  constantBodyForce_ = this->specificSettings_.template getOptionArray<double,3>("constantBodyForce", Vec3{0.0,0.0,0.0});

  if (!useAnalyticJacobian_ && !useNumericJacobian_)
  {
    LOG(WARNING) << "Cannot set both \"useAnalyticJacobian\" and \"useNumericJacobian\" to False, now using numeric jacobian.";
    useNumericJacobian_ = true;
  }

  // parse material parameters
  specificSettings_.getOptionVector("materialParameters", materialParameters_);

  LOG(DEBUG) << "HyperelasticitySolver: parsed parameters " << materialParameters_;

  if (materialParameters_.size() < Term::nMaterialParameters)
  {
    LOG(FATAL) << "Not enough material parameters specified. Specified parameters: " << materialParameters_.size()
      << " (" << materialParameters_ << "), needed parameters by " << StringUtility::demangle(typeid(Term).name()) << ": "
      << Term::nMaterialParameters;
  }

  // initialize material parameters
  // set all PARAM(i) values to the values given by materialParameters
  SEMT::set_parameters<Term::nMaterialParameters>::to(materialParameters_);

  displacementsScalingFactor_ = specificSettings_.getOptionDouble("displacementsScalingFactor", 1.0);
  dumpDenseMatlabVariables_ = specificSettings_.getOptionBool("dumpDenseMatlabVariables", false);

  // for the dynamic equation
  if (nDisplacementComponents == 6)
  {
    density_ = specificSettings_.getOptionDouble("density", 1.0, PythonUtility::Positive);
    timeStepWidth_ = specificSettings_.getOptionDouble("timeStepWidth", 1.0, PythonUtility::Positive);
  }

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_, this->specificSettings_);
  this->outputWriterManagerPressure_.initialize(this->context_["pressure"], this->context_["pressure"].getPythonConfig());
}

template<typename Term,int nDisplacementComponents>
void HyperelasticitySolver<Term,nDisplacementComponents>::
advanceTimeSpan()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  LOG(TRACE) << "advanceTimeSpan, endTime: " << endTime_;

  // write reference output values
  this->outputWriterManager_.writeOutput(this->data_, 0, 0.0);
  this->outputWriterManagerPressure_.writeOutput(this->pressureDataCopy_, 0, 0.0);

  nonlinearSolve();
  postprocessSolution();

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);

  // write current output values
  this->outputWriterManager_.writeOutput(this->data_, 1, endTime_);
  this->outputWriterManagerPressure_.writeOutput(this->pressureDataCopy_, 1, endTime_);
}

template<typename Term,int nDisplacementComponents>
void HyperelasticitySolver<Term,nDisplacementComponents>::
run()
{
  // initialize everything
  LOG(DEBUG) << "call initialize in run()";
  initialize();

  this->advanceTimeSpan();
}

template<typename Term,int nDisplacementComponents>
void HyperelasticitySolver<Term,nDisplacementComponents>::
setTimeSpan(double startTime, double endTime)
{
  endTime_ = endTime;
}

template<typename Term,int nDisplacementComponents>
void HyperelasticitySolver<Term,nDisplacementComponents>::
initialize()
{
  if (this->initialized_)
    return;

  LOG(DEBUG) << "initialize HyperelasticitySolver";
  assert(this->specificSettings_.pyObject());

  // create function space / mesh, the geometry is from the settings
  displacementsFunctionSpace_ = context_.meshManager()->functionSpace<DisplacementsFunctionSpace>(specificSettings_);

  // create 3D function space with linear basis functions
  std::vector<Vec3> nodePositionsLinearMesh;

  // loop over nodes of quadratic function space and extract nodes for linear function space
  for (node_no_t nodeIndexZ = 0; nodeIndexZ < displacementsFunctionSpace_->meshPartition()->nNodesLocalWithoutGhosts(2); nodeIndexZ++)
  {
    for (node_no_t nodeIndexY = 0; nodeIndexY < displacementsFunctionSpace_->meshPartition()->nNodesLocalWithoutGhosts(1); nodeIndexY++)
    {
      for (node_no_t nodeIndexX = 0; nodeIndexX < displacementsFunctionSpace_->meshPartition()->nNodesLocalWithoutGhosts(0); nodeIndexX++)
      {
        node_no_t nodeNoLocal = nodeIndexZ*displacementsFunctionSpace_->meshPartition()->nNodesLocalWithoutGhosts(1)*displacementsFunctionSpace_->meshPartition()->nNodesLocalWithoutGhosts(0)
          + nodeIndexY*displacementsFunctionSpace_->meshPartition()->nNodesLocalWithoutGhosts(0) + nodeIndexX;

        dof_no_t dofNoLocal = nodeNoLocal;   // no Hermite, i.e. no multiple dofs per node

        if (nodeIndexX % 2 == 0 && nodeIndexY % 2 == 0 && nodeIndexZ % 2 == 0)
        {
          nodePositionsLinearMesh.push_back(displacementsFunctionSpace_->geometryField().getValue(dofNoLocal));
        }
      }
    }
  }

  std::array<element_no_t,3> nElementsPerCoordinateDirection;
  std::array<int,3> nRanksPerCoordinateDirection;

  for (int i = 0; i < 3; i++)
  {
    nElementsPerCoordinateDirection[i] = displacementsFunctionSpace_->meshPartition()->nElementsLocal(i);
    nRanksPerCoordinateDirection[i] = displacementsFunctionSpace_->meshPartition()->nRanks(i);
  }

  pressureFunctionSpace_ = context_.meshManager()->createFunctionSpace<PressureFunctionSpace>("pressureFunctionSpace",nodePositionsLinearMesh, nElementsPerCoordinateDirection, nRanksPerCoordinateDirection);

  // initialize the data object
  // store mesh in data
  LOG(DEBUG) << "call setDisplacementsFunctionSpace of data";
  data_.setDisplacementsFunctionSpace(displacementsFunctionSpace_);
  data_.setPressureFunctionSpace(pressureFunctionSpace_);

  data_.initialize();
  pressureDataCopy_.initialize(data_.pressure(), data_.displacementsLinearMesh(), data_.velocitiesLinearMesh());
  pressureDataCopy_.setFunctionSpace(pressureFunctionSpace_);

  // create nonlinear solver PETSc context (snes)
  nonlinearSolver_ = this->context_.solverManager()->template solver<Solver::Nonlinear>(
    this->specificSettings_, this->displacementsFunctionSpace_->meshPartition()->mpiCommunicator());

  // initialize Dirichlet boundary conditions
  if (dirichletBoundaryConditions_ == nullptr)
  {
    dirichletBoundaryConditions_ = std::make_shared<DirichletBoundaryConditions<DisplacementsFunctionSpace,nDisplacementComponents>>(this->context_);
    dirichletBoundaryConditions_->initialize(this->specificSettings_, this->data_.functionSpace(), "dirichletBoundaryConditions");
  }

  // initialize Neumann boundary conditions
  if (neumannBoundaryConditions_ == nullptr)
  {
    typedef Quadrature::Gauss<3> QuadratureType;
    neumannBoundaryConditions_ = std::make_shared<NeumannBoundaryConditions<DisplacementsFunctionSpace,QuadratureType,3>>(this->context_);
    neumannBoundaryConditions_->initialize(this->specificSettings_, this->data_.functionSpace(), "neumannBoundaryConditions");
  }

  // initialize fiber direction field
  if (Term::usesFiberDirection)
  {
    initializeFiberDirections();
  }

  // setup Petsc variables
  LOG(DEBUG) << "initialize Petsc Variables";
  initializePetscVariables();

  // add this solver to the solvers diagram
  DihuContext::solverStructureVisualizer()->addSolver("HyperelasticitySolver");

  LOG(DEBUG) << "initialization done";
  this->initialized_ = true;
}

template<typename Term,int nDisplacementComponents>
void HyperelasticitySolver<Term,nDisplacementComponents>::
initializeFiberDirections()
{
  std::vector<std::string> fiberMeshNames;
  this->specificSettings_.template getOptionVector<std::string>("fiberMeshNames", fiberMeshNames);

  // loop over fiber mesh names
  for (int fiberNo = 0; fiberNo < fiberMeshNames.size(); fiberNo++)
  {
    // get fiber function space
    std::string fiberMeshName = fiberMeshNames[fiberNo];
    LOG(DEBUG) << "fiber " << fiberNo << "/" << fiberMeshNames.size() << ", mesh \"" << fiberMeshName << "\".";

    std::shared_ptr<FiberFunctionSpace> fiberFunctionSpace = context_.meshManager()->functionSpace<FiberFunctionSpace>(fiberMeshName);

    LOG(DEBUG) << "create mapping";

    // initialize mapping 1D fiber -> 3D
    context_.meshManager()->createMappingBetweenMeshes<FiberFunctionSpace,DisplacementsFunctionSpace>(
      fiberFunctionSpace, this->displacementsFunctionSpace_);
  }

  // prepare the target mesh for the mapping, set all factors to zero
  DihuContext::meshManager()->template prepareMappingLowToHigh<DisplacementsFieldVariableType>(this->data_.fiberDirection());

  // loop over fiber mesh names
  for (int fiberNo = 0; fiberNo < fiberMeshNames.size(); fiberNo++)
  {
    std::string fiberMeshName = fiberMeshNames[fiberNo];
    LOG(DEBUG) << "mesh \"" << fiberMeshName << "\".";

    std::shared_ptr<FiberFunctionSpace> fiberFunctionSpace = context_.meshManager()->functionSpace<FiberFunctionSpace>(fiberMeshName);

    // define direction field variable on the fiber that will store the direction of the fiber
    std::vector<std::string> components({"x","y","z"});
    std::shared_ptr<FieldVariable::FieldVariable<FiberFunctionSpace,3>> direction
      = fiberFunctionSpace->createFieldVariable<3>("direction", components);

    // set values of direction field variable
    int nDofsLocalWithoutGhosts = fiberFunctionSpace->nDofsLocalWithoutGhosts();
    std::vector<Vec3> geometryFieldValues;
    fiberFunctionSpace->geometryField().getValuesWithoutGhosts(geometryFieldValues);

    std::vector<Vec3> directionValues(nDofsLocalWithoutGhosts);

    // loop over local nodes
    for (dof_no_t dofNoLocal = 0; dofNoLocal < nDofsLocalWithoutGhosts; dofNoLocal++)
    {
      dof_no_t index0 = std::max((dof_no_t)(0), dofNoLocal-1);
      dof_no_t index1 = std::min((dof_no_t)(nDofsLocalWithoutGhosts-1), (dof_no_t)(dofNoLocal+1));

      // get direction of 1D fiber
      Vec3 fiberDirection = -geometryFieldValues[index0] + geometryFieldValues[index1];
      directionValues[dofNoLocal] = fiberDirection;
    }

    direction->setValuesWithoutGhosts(directionValues);


    // transfer direction values from 1D fibers to 3D field, -1 means all components
    DihuContext::meshManager()->mapLowToHighDimension<FieldVariable::FieldVariable<FiberFunctionSpace,3>, DisplacementsFieldVariableType>(
      direction, -1, this->data_.fiberDirection(), -1);
  }

  // finalize the mapping to the target mesh, compute final values by dividing by the factors
  DihuContext::meshManager()->template finalizeMappingLowToHigh<DisplacementsFieldVariableType>(this->data_.fiberDirection());

  LOG(DEBUG) << "normalize fiber direction";

  // normalize entries
  std::vector<Vec3> valuesLocalWithoutGhosts;
  this->data_.fiberDirection()->getValuesWithoutGhosts(valuesLocalWithoutGhosts);

  for (dof_no_t dofNoLocal = 0; dofNoLocal < this->data_.fiberDirection()->nDofsLocalWithoutGhosts(); dofNoLocal++)
  {
    if (MathUtility::normSquared<3>(valuesLocalWithoutGhosts[dofNoLocal]) < 1e-10)
    {
      valuesLocalWithoutGhosts[dofNoLocal] = Vec3({0, 0, 1});
      LOG(DEBUG) << "dof " << dofNoLocal << ", fiberDirection was (0,0,0), set to (0,0,1)";
    }
    else
    {
      MathUtility::normalize<3>(valuesLocalWithoutGhosts[dofNoLocal]);
      LOG(DEBUG) << "dof " << dofNoLocal << ", fiberDirection normalized: " << valuesLocalWithoutGhosts[dofNoLocal];
    }
  }

  this->data_.fiberDirection()->setValuesWithoutGhosts(valuesLocalWithoutGhosts);
}

template<typename Term,int nDisplacementComponents>
std::shared_ptr<typename HyperelasticitySolver<Term,nDisplacementComponents>::VecHyperelasticity>
HyperelasticitySolver<Term,nDisplacementComponents>::
createPartitionedPetscVec(std::string name)
{
  LOG(DEBUG) << "createPartitionedPetscVec(" << name << ")";
  return std::make_shared<VecHyperelasticity>(
    displacementsFunctionSpace_->meshPartition(), pressureFunctionSpace_->meshPartition(), dirichletBoundaryConditions_, name);
}

template<typename Term,int nDisplacementComponents>
std::shared_ptr<typename HyperelasticitySolver<Term,nDisplacementComponents>::MatHyperelasticity>
HyperelasticitySolver<Term,nDisplacementComponents>::
createPartitionedPetscMat(std::string name)
{
  // determine number of non zero entries in matrix
  int diagonalNonZeros, offdiagonalNonZeros;
  ::Data::FiniteElementsBase<DisplacementsFunctionSpace,1>::getPetscMemoryParameters(diagonalNonZeros, offdiagonalNonZeros);

  return std::make_shared<MatHyperelasticity>(
    combinedVecSolution_, 4*diagonalNonZeros, 4*offdiagonalNonZeros, name);
}


//! get the precomputed external virtual work
template<typename Term,int nDisplacementComponents>
Vec HyperelasticitySolver<Term,nDisplacementComponents>::
externalVirtualWork()
{
  return externalVirtualWorkDead_;
}

template<typename Term,int nDisplacementComponents>
void HyperelasticitySolver<Term,nDisplacementComponents>::
initializePetscVariables()
{
  /*
   * jacobian matrix layout (for one process and nDisplacementComponents==3):
   *  (U U U P)
   *  (U U U P)
   *  (U U U P)
   *  (P P P 0)
   */
  // prepare for data structures without ghost dofs, these are normal Petsc Mat's and Vec's for which all solvers are available

  // solution vector
  combinedVecSolution_ = createPartitionedPetscVec("combinedSolution");

  // residual vector
  combinedVecResidual_ = createPartitionedPetscVec("combinedResidual");

  // vector for the external virtual work contribution that does not depend on u, δW_ext,dead (this is the same as δW_ext for static case)
  combinedVecExternalVirtualWorkDead_ = createPartitionedPetscVec("combinedVecExternalVirtualWorkDead");

  // get number of local entries of the vectors, this will be the number of rows and columns of the matrix
  int nMatrixRowsLocal = combinedVecSolution_->nEntriesLocal();

  // output
  LOG(DEBUG) << "n dofs displacements: " << nDisplacementComponents << "*" << displacementsFunctionSpace_->nDofsLocalWithoutGhosts()
    << " = " << nDisplacementComponents*displacementsFunctionSpace_->nDofsLocalWithoutGhosts();
  LOG(DEBUG) << "n dofs pressure: " << pressureFunctionSpace_->nDofsLocalWithoutGhosts() << ", total: "
    << displacementsFunctionSpace_->nDofsLocalWithoutGhosts() * nDisplacementComponents + pressureFunctionSpace_->nDofsLocalWithoutGhosts();
  LOG(DEBUG) << "n BC dofs: ";

  LOG(DEBUG) << "  ux: " << dirichletBoundaryConditions_->boundaryConditionsByComponent()[0].dofNosLocal.size();
  LOG(DEBUG) << "  uy: " << dirichletBoundaryConditions_->boundaryConditionsByComponent()[1].dofNosLocal.size();
  LOG(DEBUG) << "  uz: " << dirichletBoundaryConditions_->boundaryConditionsByComponent()[2].dofNosLocal.size();

  if (nDisplacementComponents == 6)
  {
    LOG(DEBUG) << "  vx: " << dirichletBoundaryConditions_->boundaryConditionsByComponent()[3].dofNosLocal.size();
    LOG(DEBUG) << "  vy: " << dirichletBoundaryConditions_->boundaryConditionsByComponent()[4].dofNosLocal.size();
    LOG(DEBUG) << "  vz: " << dirichletBoundaryConditions_->boundaryConditionsByComponent()[5].dofNosLocal.size();
  }
  LOG(DEBUG) << "number of non-BC dofs total: " << nMatrixRowsLocal;

  // create matrix with same dof mapping as vectors
  std::shared_ptr<FunctionSpace::Generic> genericFunctionSpace = context_.meshManager()->createGenericFunctionSpace(nMatrixRowsLocal, displacementsFunctionSpace_->meshPartition(), "genericMesh");

  combinedMatrixJacobian_ = createPartitionedPetscMat("combinedJacobian");

  solverMatrixAdditionalNumericJacobian_ = PETSC_NULL;

  // if both numeric and analytic jacobian are used, create additional matrix that will hold the numeric jacobian
  if (useNumericJacobian_ && useAnalyticJacobian_)
  {
    combinedMatrixAdditionalNumericJacobian_ = createPartitionedPetscMat("combinedJacobianNumeric");

    solverMatrixAdditionalNumericJacobian_ = combinedMatrixAdditionalNumericJacobian_->valuesGlobal();
  }

  // extract the Petsc Vec's of the PartitionedPetscVecForHyperelasticity objects
  LOG(DEBUG) << "get the internal vectors";
  solverMatrixJacobian_ = combinedMatrixJacobian_->valuesGlobal();
  solverVariableSolution_ = combinedVecSolution_->valuesGlobal();
  solverVariableResidual_ = combinedVecResidual_->valuesGlobal();
  externalVirtualWorkDead_ = combinedVecExternalVirtualWorkDead_->valuesGlobal();

  // create vector with all zeros in it, this is needed for zeroing the diagonal of the stiffness matrix for initialization in evaluateAnalyticJacobian
  PetscErrorCode ierr;
  ierr = VecDuplicate(solverVariableResidual_, &zeros_); CHKERRV(ierr);
  ierr = VecZeroEntries(zeros_); CHKERRV(ierr);

  LOG(DEBUG) << "for debugging: " << combinedVecSolution_->getString();

  // compute the external virtual work, because it is constant throughout the solution process
  materialComputeExternalVirtualWorkDead();

  // output pointer values for debugging
  LOG(DEBUG) << "-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.";
  LOG(DEBUG) << "pointer value solverVariableResidual_: " << solverVariableResidual_;
  LOG(DEBUG) << "pointer value solverVariableSolution_: " << solverVariableSolution_;
  LOG(DEBUG) << "pointer value solverMatrixJacobian_:   " << solverMatrixJacobian_;
  LOG(DEBUG) << "pointer value solverMatrixAdditionalNumericJacobian_:  " << solverMatrixAdditionalNumericJacobian_;
  LOG(DEBUG) << "-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.";

  // if a numeric jacobian is in use, assemble the matrix
  // depending on if also the analytic jacobian will be computed, the storage for the numeric jacobian is either combinedMatrixJacobian_ or combinedMatrixAdditionalNumericJacobian_
  if (useNumericJacobian_)
  {
    if (useAnalyticJacobian_)
    {
      // assemble matrix, new nonzeros are allocated later
      combinedMatrixAdditionalNumericJacobian_->assembly(MAT_FINAL_ASSEMBLY);
    }
    else
    {
      // assemble matrix, new nonzeros are allocated later
      combinedMatrixJacobian_->assembly(MAT_FINAL_ASSEMBLY);
    }
  }

  if (useAnalyticJacobian_)
  {
    MatSetOption(solverMatrixJacobian_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

    // assemble matrix and nonzeros structure
    evaluateAnalyticJacobian(solverVariableSolution_, solverMatrixJacobian_);

    // output the jacobian matrix for debugging
    LOG(DEBUG) << "initial analytic jacobian matrix: ";
    dumpJacobianMatrix(solverMatrixJacobian_);
  }

  // assign all callback functions
  this->initializePetscCallbackFunctions();

  // set solution vector to zero
  this->initializeSolutionVariable();
}

//! get the PartitionedPetsVec for the residual and result of the nonlinear function
template<typename Term,int nDisplacementComponents>
std::shared_ptr<typename HyperelasticitySolver<Term,nDisplacementComponents>::VecHyperelasticity> HyperelasticitySolver<Term,nDisplacementComponents>::
combinedVecResidual()
{
  return this->combinedVecResidual_;
}

//! get the PartitionedPetsVec for the solution
template<typename Term,int nDisplacementComponents>
std::shared_ptr<typename HyperelasticitySolver<Term,nDisplacementComponents>::VecHyperelasticity> HyperelasticitySolver<Term,nDisplacementComponents>::
combinedVecSolution()
{
  return this->combinedVecSolution_;
}

template<typename Term,int nDisplacementComponents>
void HyperelasticitySolver<Term,nDisplacementComponents>::reset()
{
  this->initialized_ = false;
}

template<typename Term,int nDisplacementComponents>
typename HyperelasticitySolver<Term,nDisplacementComponents>::Data &HyperelasticitySolver<Term,nDisplacementComponents>::
data()
{
  return data_;
}

} // namespace SpatialDiscretization
