#include "specialized_solver/solid_mechanics/hyperelasticity/hyperelasticity_solver.h"

#include <Python.h>  // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "solver/solver_manager.h"
#include "data_management/specialized_solver/multidomain.h"
#include "control/diagnostic_tool/performance_measurement.h"
#include "control/diagnostic_tool/solver_structure_visualizer.h"
#include "partition/mesh_partition/01_mesh_partition_structured.h"

namespace SpatialDiscretization
{

template<typename Term,typename MeshType,int nDisplacementComponents>
HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
HyperelasticitySolver(DihuContext context, std::string settingsKey) :
  context_(context[settingsKey]), data_(context_), pressureDataCopy_(context_), initialized_(false),
  endTime_(0), lastNorm_(0), secondLastNorm_(0), currentLoadFactor_(1.0), lastSolveSucceeded_(true)
{
  // get python config
  this->specificSettings_ = this->context_.getPythonConfig();

  // read in the durationLogKey
  if (specificSettings_.hasKey("durationLogKey"))
  {
    this->durationLogKey_ = specificSettings_.getOptionString("durationLogKey", "");
  }

  // parse options concerning jacobian
  useAnalyticJacobian_  = this->specificSettings_.getOptionBool("useAnalyticJacobian", true);
  useNumericJacobian_   = this->specificSettings_.getOptionBool("useNumericJacobian", true);
  nNonlinearSolveCalls_ = this->specificSettings_.getOptionInt("nNonlinearSolveCalls", 1, PythonUtility::Positive);
  loadFactorGiveUpThreshold_ = this->specificSettings_.getOptionDouble("loadFactorGiveUpThreshold", 1e-5, PythonUtility::Positive);

  // parse constant body force, a value of "None" yields the default value, (0,0,0)
  constantBodyForce_ = this->specificSettings_.template getOptionArray<double,3>("constantBodyForce", Vec3{0.0,0.0,0.0});

  if (!useAnalyticJacobian_ && !useNumericJacobian_)
  {
    LOG(WARNING) << "Cannot set both \"useAnalyticJacobian\" and \"useNumericJacobian\" to False, now using numeric jacobian.";
    useNumericJacobian_ = true;
  }

  // parse material parameters
  specificSettings_.getOptionVector("materialParameters", materialParameters_);

  LOG(DEBUG) << "HyperelasticitySolver: parsed parameters " << materialParameters_;

  if (materialParameters_.size() != Term::nMaterialParameters)
  {
    LOG(FATAL) << "A wrong number of material parameters was specified. Specified parameters: " << materialParameters_.size()
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
    density_                 = specificSettings_.getOptionDouble("density", 1.0, PythonUtility::Positive);
    timeStepWidth_           = specificSettings_.getOptionDouble("timeStepWidth", 1.0, PythonUtility::Positive);
    extrapolateInitialGuess_ = specificSettings_.getOptionBool("extrapolateInitialGuess", true);
  }

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_, this->specificSettings_);
  this->outputWriterManagerPressure_.initialize(this->context_["pressure"], this->context_["pressure"].getPythonConfig());
  this->outputWriterManagerLoadIncrements_.initialize(this->context_["LoadIncrements"], this->context_["LoadIncrements"].getPythonConfig());
}

template<typename Term,typename MeshType,int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
advanceTimeSpan()
{
  // start duration measurement, the name of the output variable can be set by "durationLogKey" in the config
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::start(this->durationLogKey_);

  LOG(TRACE) << "advanceTimeSpan, endTime: " << endTime_;

  // write reference output values but don't increment counter
  this->outputWriterManager_.writeOutput(this->data_, 0, 0.0, 0);
  this->outputWriterManagerPressure_.writeOutput(this->pressureDataCopy_, 0, 0.0, 0);

  nonlinearSolve();
  postprocessSolution();

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);

  // write current output values
  this->outputWriterManager_.writeOutput(this->data_, 1, endTime_);
  this->outputWriterManagerPressure_.writeOutput(this->pressureDataCopy_, 1, endTime_);
}

template<typename Term,typename MeshType,int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
run()
{
  // initialize everything
  LOG(DEBUG) << "call initialize in run()";
  initialize();

  this->advanceTimeSpan();
}

template<typename Term,typename MeshType,int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
setTimeSpan(double startTime, double endTime)
{
  endTime_ = endTime;
}

template<typename Term,typename MeshType,int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
initialize()
{
  if (this->initialized_)
    return;

  LOG(DEBUG) << "initialize HyperelasticitySolver";
  assert(this->specificSettings_.pyObject());

  // create function space / mesh, the geometry is from the settings
  displacementsFunctionSpace_ = context_.meshManager()->functionSpace<DisplacementsFunctionSpace>(specificSettings_);

  pressureFunctionSpace_ = PressureFunctionSpaceCreator<typename PressureFunctionSpace::Mesh>::createPressureFunctionSpace(context_.meshManager(), displacementsFunctionSpace_);

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

  // parse load factors
  this->specificSettings_.getOptionVector("loadFactors", loadFactors_);

  if (!loadFactors_.empty())
  {
    if (fabs(loadFactors_.back() - 1.0) > 1e-12)
    {
      LOG(WARNING) << this->specificSettings_ << "[\"loadFactors\"]: Last load factor " << loadFactors_.back() << " is not 1.0.";
    }
  }

  // prepare load factors
  if (loadFactors_.empty())
    loadFactors_.push_back(1.0);

  // add this solver to the solvers diagram
  DihuContext::solverStructureVisualizer()->addSolver("HyperelasticitySolver");

  LOG(DEBUG) << "initialization done";
  this->initialized_ = true;
}

template<typename Term,typename MeshType,int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
initializeFiberDirections()
{
  std::vector<std::string> fiberMeshNames;
  this->specificSettings_.template getOptionVector<std::string>("fiberMeshNames", fiberMeshNames);

  LOG(INFO) << "initializing fiber directions...";

  LOG(DEBUG) << "config: " << this->specificSettings_;
  LOG(DEBUG) << fiberMeshNames.size() << " fiberMeshNames: " << fiberMeshNames;

  // if no fiber meshes were specified, use the settings fiberDirection
  if (fiberMeshNames.empty())
  {
    Vec3 fiberDirection = this->specificSettings_.template getOptionArray<double,3>("fiberDirection", Vec3{0,0,1});

    LOG(DEBUG) << "displacements field variable type data: " << StringUtility::demangle(typeid(typename Data::DisplacementsFieldVariableType).name());
    LOG(DEBUG) << "displacements function space type: " << StringUtility::demangle(typeid(DisplacementsFunctionSpace).name());
    LOG(DEBUG) << "displacements function space mesh partition: " << *this->displacementsFunctionSpace_->meshPartition();

    std::vector<Vec3> fiberDirections(this->displacementsFunctionSpace_->nDofsLocalWithGhosts(), fiberDirection);
    this->data_.fiberDirection()->setValues(this->displacementsFunctionSpace_->meshPartition()->dofNosLocal(), fiberDirections);
  }
  else
  {
    // loop over fiber mesh names
    for (int fiberNo = 0; fiberNo < fiberMeshNames.size(); fiberNo++)
    {
      // get fiber function space
      std::string fiberMeshName = fiberMeshNames[fiberNo];
      LOG(DEBUG) << "fiber " << fiberNo << "/" << fiberMeshNames.size() << ", mesh \"" << fiberMeshName << "\".";

      std::shared_ptr<FiberFunctionSpace> fiberFunctionSpace = context_.meshManager()->functionSpace<FiberFunctionSpace>(fiberMeshName);

      LOG(DEBUG) << "create mapping  1D fiber -> 3D: \"" << fiberFunctionSpace->meshName() << "\" -> \"" << this->displacementsFunctionSpace_->meshName() << "\".";

      // initialize mapping 1D fiber -> 3D
      DihuContext::mappingBetweenMeshesManager()->template mappingBetweenMeshes<FiberFunctionSpace,DisplacementsFunctionSpace>(fiberFunctionSpace, this->displacementsFunctionSpace_);
    }

    using SourceFunctionSpaceType = FieldVariable::FieldVariable<FiberFunctionSpace,3>;
    using TargetFunctionSpaceType = DisplacementsFieldVariableType;

    // prepare the target mesh for the mapping, set all factors to zero
    std::shared_ptr<FieldVariable::FieldVariable<FiberFunctionSpace,3>> direction;
    DihuContext::mappingBetweenMeshesManager()->template prepareMapping<SourceFunctionSpaceType,TargetFunctionSpaceType>(direction, this->data_.fiberDirection(), -1);

    // loop over fiber mesh names
    for (int fiberNo = 0; fiberNo < fiberMeshNames.size(); fiberNo++)
    {
      std::string fiberMeshName = fiberMeshNames[fiberNo];
      LOG(DEBUG) << "mesh \"" << fiberMeshName << "\".";

      std::shared_ptr<FiberFunctionSpace> fiberFunctionSpace = context_.meshManager()->functionSpace<FiberFunctionSpace>(fiberMeshName);

      // define direction field variable on the fiber that will store the direction of the fiber
      std::vector<std::string> components({"x","y","z"});
      direction = fiberFunctionSpace->createFieldVariable<3>("direction", components);

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
      DihuContext::mappingBetweenMeshesManager()->template map<SourceFunctionSpaceType,TargetFunctionSpaceType>(
        direction, this->data_.fiberDirection(), -1, -1, false);
    }

    // finalize the mapping to the target mesh, compute final values by dividing by the factors
    DihuContext::mappingBetweenMeshesManager()->template finalizeMapping<SourceFunctionSpaceType,TargetFunctionSpaceType>(direction, this->data_.fiberDirection(), -1, -1, false);
  }

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

  this->data_.fiberDirection()->zeroGhostBuffer();
  this->data_.fiberDirection()->finishGhostManipulation();
}

template<typename Term,typename MeshType,int nDisplacementComponents>
std::shared_ptr<typename HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::VecHyperelasticity>
HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
createPartitionedPetscVec(std::string name)
{
  LOG(DEBUG) << "createPartitionedPetscVec(" << name << ")";
  return std::make_shared<VecHyperelasticity>(
    displacementsFunctionSpace_->meshPartition(), pressureFunctionSpace_->meshPartition(), dirichletBoundaryConditions_, name);
}

template<typename Term,typename MeshType,int nDisplacementComponents>
std::shared_ptr<typename HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::MatHyperelasticity>
HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
createPartitionedPetscMat(std::string name)
{
  // determine number of non zero entries in matrix
  int nNonZerosDiagonal, nNonZerosOffdiagonal;
  ::Data::FiniteElementsBase<DisplacementsFunctionSpace,1>::getPetscMemoryParameters(nNonZerosDiagonal, nNonZerosOffdiagonal);

  return std::make_shared<MatHyperelasticity>(
    combinedVecSolution_, 4*nNonZerosDiagonal, 4*nNonZerosOffdiagonal, name);
}


//! get the precomputed external virtual work
template<typename Term,typename MeshType,int nDisplacementComponents>
Vec HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
externalVirtualWork()
{
  return externalVirtualWorkDead_;
}

template<typename Term,typename MeshType,int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
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
  ierr = VecDuplicate(solverVariableResidual_, &lastSolution_); CHKERRV(ierr);

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
template<typename Term,typename MeshType,int nDisplacementComponents>
std::shared_ptr<typename HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::VecHyperelasticity> HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
combinedVecResidual()
{
  return this->combinedVecResidual_;
}

//! get the PartitionedPetsVec for the solution
template<typename Term,typename MeshType,int nDisplacementComponents>
std::shared_ptr<typename HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::VecHyperelasticity> HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
combinedVecSolution()
{
  return this->combinedVecSolution_;
}

template<typename Term,typename MeshType,int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::reset()
{
  this->initialized_ = false;
}

template<typename Term,typename MeshType,int nDisplacementComponents>
typename HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::Data &HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
data()
{
  return data_;
}

//! get a pointer to the dirichlet boundary conditions object
template<typename Term,typename MeshType,int nDisplacementComponents>
std::shared_ptr<DirichletBoundaryConditions<typename HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::DisplacementsFunctionSpace,nDisplacementComponents>>
HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
dirichletBoundaryConditions()
{
  return dirichletBoundaryConditions_;
}

//! set new neumann bc's = traction for the next solve
template<typename Term,typename MeshType,int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
updateNeumannBoundaryConditions(std::shared_ptr<NeumannBoundaryConditions<typename HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::DisplacementsFunctionSpace,Quadrature::Gauss<3>,3>> newNeumannBoundaryConditions)
{
  neumannBoundaryConditions_ = newNeumannBoundaryConditions;

  // compute new value for the rhs, δW_ext,dead = int_Ω B^L * phi^L * phi^M * δu^M dx + int_∂Ω T^L * phi^L * phi^M * δu^M dS
  materialComputeExternalVirtualWorkDead();
}

template<typename Term,typename MeshType,int nDisplacementComponents>
void HyperelasticitySolver<Term,MeshType,nDisplacementComponents>::
addDirichletBoundaryConditions(std::vector<typename DirichletBoundaryConditions<DisplacementsFunctionSpace,nDisplacementComponents>::ElementWithNodes> &boundaryConditionElements, bool overwriteBcOnSameDof)
{
  dirichletBoundaryConditions_->addBoundaryConditions(boundaryConditionElements, overwriteBcOnSameDof);

  // save previous values
  std::array<std::vector<double>, nDisplacementComponents+1> combinedSolutionValues;
  std::array<std::vector<double>, nDisplacementComponents+1> combinedResidualValues;
  std::array<std::vector<double>, nDisplacementComponents+1> combinedVecExternalVirtualWorkDeadValues;

  int nDisplacementsVelocityValues = this->displacementsFunctionSpace_->meshPartition()->nDofsLocalWithoutGhosts();
  int nPressureValues = this->pressureFunctionSpace_->meshPartition()->nDofsLocalWithoutGhosts();
  std::vector<int> indices(nDisplacementsVelocityValues);
  std::iota(indices.begin(), indices.end(), 0);

  // loop over displacement components and one pressure component
  for (int componentNo = 0; componentNo < nDisplacementComponents+1; componentNo++)
  {
    int nValues = nDisplacementsVelocityValues;
    if (componentNo == nDisplacementComponents)
      nValues = nPressureValues;

    // solution vector
    combinedSolutionValues[componentNo].resize(nValues);
    combinedVecSolution_->getValues(componentNo, nValues, indices.data(), combinedSolutionValues[componentNo].data());

    // residual vector
    combinedResidualValues[componentNo].resize(nValues);
    combinedVecResidual_->getValues(componentNo, nValues, indices.data(), combinedResidualValues[componentNo].data());

    // vector for the external virtual work contribution that does not depend on u, δW_ext,dead (this is the same as δW_ext for static case)
    combinedVecExternalVirtualWorkDeadValues[componentNo].resize(nValues);
    combinedVecExternalVirtualWorkDead_->getValues(componentNo, nValues, indices.data(), combinedVecExternalVirtualWorkDeadValues[componentNo].data());
  }

  // remove generic meshes "genericMesh" and "genericMeshForMatrixcombinedJacobian"
  DihuContext::meshManager()->deleteFunctionSpace("genericMesh");
  DihuContext::meshManager()->deleteFunctionSpace("genericMeshForMatrixcombinedJacobian");

  // create new vectors and matrices with the updated boundary conditions
  initializePetscVariables();

  // restore previous values, except for the new Dirichlet values
  // loop over displacement components and one pressure component
  for (int componentNo = 0; componentNo < nDisplacementComponents+1; componentNo++)
  {
    int nValues = nDisplacementsVelocityValues;
    if (componentNo == nDisplacementComponents)
      nValues = nPressureValues;

    // solution vector
    combinedSolutionValues[componentNo].resize(nValues);
    combinedVecSolution_->setValues(componentNo, nValues, indices.data(), combinedSolutionValues[componentNo].data());

    // residual vector
    combinedResidualValues[componentNo].resize(nValues);
    combinedVecResidual_->getValues(componentNo, nValues, indices.data(), combinedResidualValues[componentNo].data());

    // vector for the external virtual work contribution that does not depend on u, δW_ext,dead (this is the same as δW_ext for static case)
    combinedVecExternalVirtualWorkDeadValues[componentNo].resize(nValues);
    combinedVecExternalVirtualWorkDead_->getValues(componentNo, nValues, indices.data(), combinedVecExternalVirtualWorkDeadValues[componentNo].data());
  }

  combinedVecSolution_->startGhostManipulation();
  combinedVecSolution_->finishGhostManipulation();

  combinedVecResidual_->startGhostManipulation();
  combinedVecResidual_->finishGhostManipulation();

  PetscErrorCode ierr;
  ierr = VecAssemblyBegin(solverVariableSolution_); CHKERRV(ierr);
  ierr = VecAssemblyEnd(solverVariableSolution_); CHKERRV(ierr);
}


} // namespace SpatialDiscretization
