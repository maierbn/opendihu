#include "specialized_solver/solid_mechanics/hyperelasticity/hyperelasticity_solver.h"

#include <Python.h>  // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "data_management/specialized_solver/multidomain.h"
#include "control/performance_measurement.h"

namespace SpatialDiscretization
{

template<typename Term>
HyperelasticitySolver<Term>::
HyperelasticitySolver(DihuContext context) :
  context_(context["HyperelasticitySolver"]), data_(context_), pressureDataCopy_(context_), initialized_(false), endTime_(0)
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

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_, this->specificSettings_);
  this->outputWriterManagerPressure_.initialize(this->context_, this->context_["pressure"].getPythonConfig());
}

template<typename Term>
void HyperelasticitySolver<Term>::
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

  // stop duration measurement
  if (this->durationLogKey_ != "")
    Control::PerformanceMeasurement::stop(this->durationLogKey_);

  // write current output values
  this->outputWriterManager_.writeOutput(this->data_, 1, endTime_);
  this->outputWriterManagerPressure_.writeOutput(this->pressureDataCopy_, 1, endTime_);
}

template<typename Term>
void HyperelasticitySolver<Term>::
run()
{
  // initialize everything
  initialize();

  this->advanceTimeSpan();
}

template<typename Term>
void HyperelasticitySolver<Term>::
setTimeSpan(double startTime, double endTime)
{
  endTime_ = endTime;
}

template<typename Term>
void HyperelasticitySolver<Term>::
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
  pressureDataCopy_.initialize(data_.pressure(), data_.displacementsLinearMesh());
  pressureDataCopy_.setFunctionSpace(pressureFunctionSpace_);

  // initialize Dirichlet boundary conditions
  if (dirichletBoundaryConditions_ == nullptr)
  {
    dirichletBoundaryConditions_ = std::make_shared<DirichletBoundaryConditions<DisplacementsFunctionSpace,3>>(this->context_);
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
  initializePetscVariables();

  LOG(DEBUG) << "initialization done";
  this->initialized_ = true;
}

template<typename Term>
void HyperelasticitySolver<Term>::
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
  DihuContext::meshManager()->template prepareMapping<DisplacementsFieldVariableType>(this->data_.fiberDirection());


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
      dof_no_t index0 = std::max(0, dofNoLocal-1);
      dof_no_t index1 = std::min(nDofsLocalWithoutGhosts-1, dofNoLocal+1);

      Vec3 direction = -geometryFieldValues[index0] + geometryFieldValues[index1];
      directionValues[dofNoLocal] = direction;
    }

    direction->setValuesWithoutGhosts(directionValues);

    // transfer direction values
    DihuContext::meshManager()->mapLowToHighDimension<FieldVariable::FieldVariable<FiberFunctionSpace,3>, DisplacementsFieldVariableType>(
      direction, this->data_.fiberDirection());
  }

  // finalize the mapping to the target mesh, compute final values by dividing by the factors
  DihuContext::meshManager()->template finalizeMapping<DisplacementsFieldVariableType>(this->data_.fiberDirection());

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

template<typename Term>
void HyperelasticitySolver<Term>::
initializePetscVariables()
{
  /*
   * jacobian matrix layout (for one process):
   *  (U U U P)
   *  (U U U P)
   *  (U U U P)
   *  (P P P 0)
   */

  // determine number of non zero entries in matrix
  int diagonalNonZeros, offdiagonalNonZeros;
  ::Data::FiniteElementsBase<DisplacementsFunctionSpace,1>::getPetscMemoryParameters(diagonalNonZeros, offdiagonalNonZeros);

  // prepare for data structures without ghost dofs, these are normal Petsc Mat's and Vec's for which all solvers are available

  // solution vector
  combinedVecSolution_ = std::make_shared<PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace>>(
    displacementsFunctionSpace_->meshPartition(), pressureFunctionSpace_->meshPartition(), dirichletBoundaryConditions_, "combinedSolution");

  // residual vector
  combinedVecResidual_ = std::make_shared<PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace>>(
    displacementsFunctionSpace_->meshPartition(), pressureFunctionSpace_->meshPartition(), dirichletBoundaryConditions_, "combinedResidual");

  // vector for the external virtual work contribution, δW_ext
  combinedVecExternalVirtualWork_ = std::make_shared<PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace>>(
    displacementsFunctionSpace_->meshPartition(), pressureFunctionSpace_->meshPartition(), dirichletBoundaryConditions_, "combinedVecExternalVirtualWork");

  // get number of local entries of the vectors, this will be the number of rows and columns of the matrix
  int nMatrixRowsLocal = combinedVecSolution_->nEntriesLocal();

  // output
  LOG(DEBUG) << "n dofs displacements: 3*" << displacementsFunctionSpace_->nDofsLocalWithoutGhosts() << " = " << 3*displacementsFunctionSpace_->nDofsLocalWithoutGhosts();
  LOG(DEBUG) << "n dofs pressure: " << pressureFunctionSpace_->nDofsLocalWithoutGhosts() << ", total: " << displacementsFunctionSpace_->nDofsLocalWithoutGhosts() * 3 + pressureFunctionSpace_->nDofsLocalWithoutGhosts();
  LOG(DEBUG) << "n BC dofs: ux:" << dirichletBoundaryConditions_->boundaryConditionsByComponent()[0].dofNosLocal.size()
    << ", uy: " << dirichletBoundaryConditions_->boundaryConditionsByComponent()[1].dofNosLocal.size()
    << ", uz: " << dirichletBoundaryConditions_->boundaryConditionsByComponent()[2].dofNosLocal.size();
  LOG(DEBUG) << "number of non-BC dofs total: " << nMatrixRowsLocal;

  // create matrix with same dof mapping as vectors
  std::shared_ptr<FunctionSpace::Generic> genericFunctionSpace = context_.meshManager()->createGenericFunctionSpace(nMatrixRowsLocal, displacementsFunctionSpace_->meshPartition(), "genericMesh");

  combinedMatrixJacobian_ = std::make_shared<PartitionedPetscMatForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace>>(
    combinedVecSolution_, 4*diagonalNonZeros, 4*offdiagonalNonZeros, "combinedJacobian");

  solverMatrixAdditionalNumericJacobian_ = PETSC_NULL;

  // if both numeric and analytic jacobian are used, create additional matrix that will hold the numeric jacobian
  if (useNumericJacobian_ && useAnalyticJacobian_)
  {
    combinedMatrixAdditionalNumericJacobian_ = std::make_shared<PartitionedPetscMatForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace>>(
      combinedVecSolution_, 4*diagonalNonZeros, 4*offdiagonalNonZeros, "combinedJacobianNumeric");

    solverMatrixAdditionalNumericJacobian_ = combinedMatrixAdditionalNumericJacobian_->valuesGlobal();
  }

  // extract the Petsc Vec's of the PartitionedPetscVecForHyperelasticity objects
  LOG(DEBUG) << "get the internal vectors";
  solverMatrixJacobian_ = combinedMatrixJacobian_->valuesGlobal();
  solverVariableSolution_ = combinedVecSolution_->valuesGlobal();
  solverVariableResidual_ = combinedVecResidual_->valuesGlobal();
  externalVirtualWork_ = combinedVecExternalVirtualWork_->valuesGlobal();

  // create vector with all zeros in it, this is needed for zeroing the diagonal of the stiffness matrix for initialization in evaluateAnalyticJacobian
  PetscErrorCode ierr;
  ierr = VecDuplicate(solverVariableResidual_, &zeros_); CHKERRV(ierr);
  ierr = VecZeroEntries(zeros_); CHKERRV(ierr);

  LOG(DEBUG) << "for debugging: " << combinedVecSolution_->getString();

  // compute the external virtual work, because it is constant throughout the solution process
  combinedVecExternalVirtualWork_->zeroEntries();
  combinedVecExternalVirtualWork_->startGhostManipulation();

  // get externalVirtualWork_
  std::vector<double> values;
  for (int componentNo = 0; componentNo < 3; componentNo++)
  {
    values.clear();
    neumannBoundaryConditions_->rhs()->getValuesWithoutGhosts(componentNo, values);
    LOG(DEBUG) << "component " << componentNo << ", neumannBoundaryConditions_ rhs values: " << values;


    combinedVecExternalVirtualWork_->setValues(componentNo, displacementsFunctionSpace_->meshPartition()->nDofsLocalWithoutGhosts(),
                                                displacementsFunctionSpace_->meshPartition()->dofNosLocal().data(), values.data(), INSERT_VALUES);
  }

  combinedVecExternalVirtualWork_->finishGhostManipulation();

  LOG(DEBUG) << "combinedVecExternalVirtualWork: " << combinedVecExternalVirtualWork_->getString();
  combinedVecExternalVirtualWork_->startGhostManipulation();

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
}


//! get the PartitionedPetsVec for the residual and result of the nonlinear function
template<typename Term>
std::shared_ptr<PartitionedPetscVecForHyperelasticity<typename HyperelasticitySolver<Term>::DisplacementsFunctionSpace,typename HyperelasticitySolver<Term>::PressureFunctionSpace>> HyperelasticitySolver<Term>::
combinedVecResidual()
{
  return this->combinedVecResidual_;
}

//! get the PartitionedPetsVec for the solution
template<typename Term>
std::shared_ptr<PartitionedPetscVecForHyperelasticity<typename HyperelasticitySolver<Term>::DisplacementsFunctionSpace,typename HyperelasticitySolver<Term>::PressureFunctionSpace>> HyperelasticitySolver<Term>::
combinedVecSolution()
{
  return this->combinedVecSolution_;
}

template<typename Term>
void HyperelasticitySolver<Term>::reset()
{
  this->initialized_ = false;
}

template<typename Term>
typename HyperelasticitySolver<Term>::Data &HyperelasticitySolver<Term>::
data()
{
  return data_;
}

} // namespace SpatialDiscretization
