#include "specialized_solver/hyperelasticity/quasi_static_hyperelasticity_solver.h"

#include <Python.h>  // has to be the first included header

#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include "data_management/specialized_solver/multidomain.h"
#include "control/performance_measurement.h"

namespace TimeSteppingScheme
{

QuasiStaticHyperelasticitySolver::
QuasiStaticHyperelasticitySolver(DihuContext context) :
  context_(context["QuasiStaticHyperelasticitySolver"]), data_(context_), pressureDataCopy_(context_), initialized_(false), endTime_(0)
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
  c1_ = specificSettings_.getOptionDouble("c1", 1.0, PythonUtility::ValidityCriterion::NonNegative);
  c2_ = specificSettings_.getOptionDouble("c2", 1.0, PythonUtility::ValidityCriterion::NonNegative);
  displacementsScalingFactor_ = specificSettings_.getOptionDouble("scalingFactor", 1.0);

  // initialize material parameters
  std::vector<double> parametersVector({c1_, c2_});

  // set all PARAM(i) values to the values given by materialParameters
  SEMT::set_parameters<2>::to(parametersVector);

  LOG(DEBUG) << "QuasiStaticHyperelasticitySolver: parsed parameters c1: " << c1_ << ", c2: " << c2_;
  LOG(DEBUG) << "now parse output writers";

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_, this->specificSettings_);
  this->outputWriterManagerPressure_.initialize(this->context_, this->context_["pressure"].getPythonConfig());
}


void QuasiStaticHyperelasticitySolver::
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

void QuasiStaticHyperelasticitySolver::
run()
{
  // initialize everything
  initialize();

  this->advanceTimeSpan();
}


void QuasiStaticHyperelasticitySolver::
setTimeSpan(double startTime, double endTime)
{
  endTime_ = endTime;
}


void QuasiStaticHyperelasticitySolver::
initialize()
{
  if (this->initialized_)
    return;

  LOG(DEBUG) << "initialize QuasiStaticHyperelasticitySolver";
  assert(this->specificSettings_.pyObject());

  // create function space / mesh, the geometry is from the settings
  displacementsFunctionSpace_ = context_.meshManager()->functionSpace<DisplacementsFunctionSpace>(specificSettings_);

  // create 3D function space with linear basis functions
  //FunctionSpaceDofsNodes(std::shared_ptr<Partition::Manager> partitionManager, const std::vector<Vec3> &nodePositions,
  //                       const std::array<element_no_t,D> nElementsPerCoordinateDirection, const std::array<int,D> nRanksPerCoordinateDirection);

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
    dirichletBoundaryConditions_ = std::make_shared<SpatialDiscretization::DirichletBoundaryConditions<DisplacementsFunctionSpace,3>>(this->context_);
    dirichletBoundaryConditions_->initialize(this->specificSettings_, this->data_.functionSpace(), "dirichletBoundaryConditions");
  }

  // initialize Neumann boundary conditions
  if (neumannBoundaryConditions_ == nullptr)
  {
    typedef Quadrature::Gauss<3> QuadratureType;
    neumannBoundaryConditions_ = std::make_shared<SpatialDiscretization::NeumannBoundaryConditions<DisplacementsFunctionSpace,QuadratureType,3>>(this->context_);
    neumannBoundaryConditions_->initialize(this->specificSettings_, this->data_.functionSpace(), "neumannBoundaryConditions");
  }

  // setup Petsc variables

  // create nested matrix for the jacobian matrix for the Newton solver, which in case of nonlinear elasticity is the tangent stiffness matrix
  PetscErrorCode ierr;
  MPI_Comm mpiCommunicator = displacementsFunctionSpace_->meshPartition()->mpiCommunicator();

  /*
   *  U U U P
   *  U U U P
   *  U U U P
   *  P P P P
   */

  // determine number of non zero entries in matrix
  int diagonalNonZeros, offdiagonalNonZeros;
  ::Data::FiniteElementsBase<DisplacementsFunctionSpace,1>::getPetscMemoryParameters(diagonalNonZeros, offdiagonalNonZeros);

  // use Petsc nested matrix data structures, these avoid to copy values from the geometry fields, but only gmres and cg solvers are implemented
  if (useNestedMat_)
  {
    // top u matrices
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        std::stringstream name;
        name << "jacobian" << i << j;
        uMatrix_.emplace_back(displacementsFunctionSpace_->meshPartition(), 1, diagonalNonZeros, offdiagonalNonZeros, name.str());
        submatrices_[i*4 + j] = uMatrix_.back().valuesGlobal();
      }
    }

    // right p matrices
    for (int i = 0; i < 3; i++)
    {
      std::stringstream name;
      name << "jacobian" << i << 3;
      upMatrix_.emplace_back(displacementsFunctionSpace_->meshPartition(), pressureFunctionSpace_->meshPartition(),
                            1, diagonalNonZeros, offdiagonalNonZeros, name.str());
      submatrices_[i*4 + 3] = upMatrix_.back().valuesGlobal();
    }

    // bottom p matrices
    for (int j = 0; j < 3; j++)
    {
      std::stringstream name;
      name << "jacobian" << 3 << j;
      puMatrix_.emplace_back(pressureFunctionSpace_->meshPartition(), displacementsFunctionSpace_->meshPartition(),
                            1, diagonalNonZeros, offdiagonalNonZeros, name.str());
      submatrices_[3*4 + j] = puMatrix_.back().valuesGlobal();
    }

    // bottom right matrix
    pMatrix_.emplace_back(pressureFunctionSpace_->meshPartition(), 1, diagonalNonZeros, offdiagonalNonZeros, "jacobian33");
    submatrices_[3*4 + 3] = pMatrix_.back().valuesGlobal();

    // create 4x4 nested matrix
    ierr = MatCreateNest(mpiCommunicator, 4, NULL, 4, NULL, submatrices_.data(), &this->solverMatrixJacobian_); CHKERRV(ierr);

    // create nested vector for solution
    for (int i = 0; i < 3; i++)
    {
      subvectorsSolution_[i] = this->data_.displacements()->valuesGlobal(i);

      VecDuplicate(subvectorsSolution_[i], &subvectorsResidual_[i]);
    }
    subvectorsSolution_[3] = this->data_.pressure()->valuesGlobal();
    VecDuplicate(subvectorsSolution_[3], &subvectorsResidual_[3]);

    ierr = VecCreateNest(mpiCommunicator, 4, NULL, subvectorsSolution_.data(), &this->solverVariableSolution_); CHKERRV(ierr);

    // create nested vector for residual
    ierr = VecCreateNest(mpiCommunicator, 4, NULL, subvectorsResidual_.data(), &this->solverVariableResidual_); CHKERRV(ierr);
  }
  else
  {
    // prepare for non-nested data structures without ghost dofs, these are normal Petsc Mat's and Vec's for which all solvers are available

    combinedVecSolution_ = std::make_shared<PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace>>(
      displacementsFunctionSpace_->meshPartition(), pressureFunctionSpace_->meshPartition(), dirichletBoundaryConditions_, "combinedSolution");

    combinedVecResidual_ = std::make_shared<PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace>>(
      displacementsFunctionSpace_->meshPartition(), pressureFunctionSpace_->meshPartition(), dirichletBoundaryConditions_, "combinedResidual");

    combinedVecExternalVirtualWork_ = std::make_shared<PartitionedPetscVecForHyperelasticity<DisplacementsFunctionSpace,PressureFunctionSpace>>(
      displacementsFunctionSpace_->meshPartition(), pressureFunctionSpace_->meshPartition(), dirichletBoundaryConditions_, "combinedVecExternalVirtualWork");

    int nMatrixRowsLocal = combinedVecSolution_->nEntriesLocal();

    LOG(DEBUG) << "n dofs displacements: 3*" << displacementsFunctionSpace_->nDofsLocalWithoutGhosts() << " = " << 3*displacementsFunctionSpace_->nDofsLocalWithoutGhosts();
    LOG(DEBUG) << "n dofs pressure: " << pressureFunctionSpace_->nDofsLocalWithoutGhosts() << ", total: " << displacementsFunctionSpace_->nDofsLocalWithoutGhosts() * 3 + pressureFunctionSpace_->nDofsLocalWithoutGhosts();
    LOG(DEBUG) << "n BC dofs: ux:" << dirichletBoundaryConditions_->boundaryConditionsByComponent()[0].dofNosLocal.size()
      << ", uy: " << dirichletBoundaryConditions_->boundaryConditionsByComponent()[1].dofNosLocal.size()
      << ", uz: " << dirichletBoundaryConditions_->boundaryConditionsByComponent()[2].dofNosLocal.size();
    LOG(DEBUG) << "number of non-BC dofs total: " << nMatrixRowsLocal;

    std::shared_ptr<FunctionSpace::Generic> genericFunctionSpace = context_.meshManager()->createGenericFunctionSpace(nMatrixRowsLocal, displacementsFunctionSpace_->meshPartition(), "genericMesh");

    combinedMatrixJacobian_ = std::make_shared<PartitionedPetscMat<FunctionSpace::Generic>>(
      genericFunctionSpace->meshPartition(), 1, 4*diagonalNonZeros, 4*offdiagonalNonZeros, "combinedJacobian");

    solverMatrixJacobian_ = combinedMatrixJacobian_->valuesGlobal();
    solverVariableSolution_ = combinedVecSolution_->valuesGlobal();
    solverVariableResidual_ = combinedVecResidual_->valuesGlobal();
    externalVirtualWork_ = combinedVecExternalVirtualWork_->valuesGlobal();

    LOG(DEBUG) << "for debugging: " << combinedVecSolution_->getString();

    combinedVecExternalVirtualWork_->zeroEntries();
    combinedMatrixJacobian_->assembly(MAT_FINAL_ASSEMBLY);

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
  }

  LOG(DEBUG) << "-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.";
  LOG(DEBUG) << "pointer value solverVariableResidual_: " << solverVariableResidual_;
  LOG(DEBUG) << "pointer value solverVariableSolution_: " << solverVariableSolution_;
  LOG(DEBUG) << "pointer value solverMatrixJacobian_:   " << solverMatrixJacobian_;

  if (useNestedMat_)
  {
    for (int i = 0; i < 3; i++)
    {
      LOG(DEBUG) << "subvectorsResidual_[" << i << "]: " << subvectorsResidual_[i];
    }
  }
  LOG(DEBUG) << "-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.";

  // assemble matrix
  evaluateAnalyticJacobian(solverVariableSolution_, solverMatrixJacobian_);

  LOG(DEBUG) << "initialization done";
  this->initialized_ = true;
}

//! get the PartitionedPetsVec for the residual and result of the nonlinear function
std::shared_ptr<PartitionedPetscVecForHyperelasticity<typename QuasiStaticHyperelasticitySolver::DisplacementsFunctionSpace,typename QuasiStaticHyperelasticitySolver::PressureFunctionSpace>> QuasiStaticHyperelasticitySolver::
combinedVecResidual()
{
  return this->combinedVecResidual_;
}

//! get the PartitionedPetsVec for the solution
std::shared_ptr<PartitionedPetscVecForHyperelasticity<typename QuasiStaticHyperelasticitySolver::DisplacementsFunctionSpace,typename QuasiStaticHyperelasticitySolver::PressureFunctionSpace>> QuasiStaticHyperelasticitySolver::
combinedVecSolution()
{
  return this->combinedVecSolution_;
}

void QuasiStaticHyperelasticitySolver::reset()
{
  this->initialized_ = false;
}

//! return whether the underlying discretizableInTime object has a specified mesh type and is not independent of the mesh type

bool QuasiStaticHyperelasticitySolver::
knowsMeshType()
{
  return true;
}

typename QuasiStaticHyperelasticitySolver::Data &
QuasiStaticHyperelasticitySolver::
data()
{
  return data_;
}

} // namespace TimeSteppingScheme
