#include "control/precice/volume_coupling/muscle_contraction.h"

#include <sstream>

namespace PreciceAdapter
{

template<typename NestedSolver>
MuscleContraction<NestedSolver>::
MuscleContraction(DihuContext context) :
  Runnable(),
  context_(context["MuscleContraction"]), nestedSolver_(this->context_), initialized_(false)
{
  // get python settings object from context
  this->specificSettings_ = this->context_.getPythonConfig();
}

template<typename NestedSolver>
void MuscleContraction<NestedSolver>::
initialize()
{
#ifdef HAVE_PRECICE

  // make sure that we initialize only once, in the next call, initialized_ is true
  if (initialized_)
    return;

  // initialize() will be called before the simulation starts.

  // add this solver to the solvers diagram, which is an ASCII art representation that will be created at the end of the simulation.
  DihuContext::solverStructureVisualizer()->addSolver("PreciceAdapter::MuscleContraction", true);   // hasInternalConnectionToFirstNestedSolver=true (the last argument) means slot connector data is shared with the first subsolver

  // indicate in solverStructureVisualizer that now a child solver will be initialized
  DihuContext::solverStructureVisualizer()->beginChild();

  // call initialize of the nested solver
  nestedSolver_.initialize();

  // indicate in solverStructureVisualizer that the child solver initialization is done
  DihuContext::solverStructureVisualizer()->endChild();

  // set the slotConnectorData for the solverStructureVisualizer to appear in the solver diagram
  DihuContext::solverStructureVisualizer()->setSlotConnectorData(getSlotConnectorData());

  connectorSlotIdGamma_ = this->specificSettings_.getOptionInt("connectorSlotIdGamma", 2, PythonUtility::Positive);

  // initialize precice
  const std::string solverName = "MuscleContraction";
  const std::string configFileName = this->specificSettings_.getOptionString("preciceConfigFilename", "../precice-config.xml");

  int rankNo = nestedSolver_.data().functionSpace()->meshPartition()->rankSubset()->ownRankNo();
  int nRanks = nestedSolver_.data().functionSpace()->meshPartition()->rankSubset()->size();

  preciceSolverInterface_ = std::make_unique<precice::SolverInterface>(solverName, configFileName, rankNo, nRanks);

  // store the node positions to precice
  std::string meshName = "MuscleContractionMesh";
  preciceMeshId_ = preciceSolverInterface_->getMeshID(meshName);

  std::shared_ptr<::FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, ::BasisFunction::LagrangeOfOrder<2>>> functionSpace
    = nestedSolver_.data().functionSpace();

  int nDofsLocalWithoutGhosts = functionSpace->nDofsLocalWithoutGhosts();
  std::vector<Vec3> geometryValues;
  functionSpace->geometryField().getValuesWithoutGhosts(geometryValues);

  std::vector<double> geometryValuesPrecice(3*nDofsLocalWithoutGhosts);
  preciceVertexIds_.resize(nDofsLocalWithoutGhosts);

  for (dof_no_t dofNoLocal = 0; dofNoLocal < nDofsLocalWithoutGhosts; dofNoLocal++)
  {
    for (int i = 0; i < 3; i++)
    {
      geometryValuesPrecice[3*dofNoLocal + i] = geometryValues[dofNoLocal][i];
    }
  }

  LOG(DEBUG) << "setMeshVertices to precice, " << 3*nDofsLocalWithoutGhosts << " values";
  preciceSolverInterface_->setMeshVertices(preciceMeshId_, nDofsLocalWithoutGhosts, geometryValuesPrecice.data(), preciceVertexIds_.data());

  LOG(DEBUG) << "precice defined vertexIds: " << preciceVertexIds_;

  // initialize data ids
  preciceDataIdGeometry_  = preciceSolverInterface_->getDataID("Geometry",  preciceMeshId_);
  preciceDataIdGamma_     = preciceSolverInterface_->getDataID("Gamma",     preciceMeshId_);
  //preciceDataIdLambda_    = preciceSolverInterface_->getDataID("Lambda",    preciceMeshId_);
  //preciceDataIdLambdaDot_ = preciceSolverInterface_->getDataID("LambdaDot", preciceMeshId_);

  LOG(DEBUG) << "data id geometry: " << preciceDataIdGeometry_;
  LOG(DEBUG) << "data id gamma: " << preciceDataIdGamma_;

  //preciceSolverInterface_->setMeshVertices(preciceMeshId_, int size, double* positions, int* ids);

  maximumPreciceTimestepSize_ = preciceSolverInterface_->initialize();

  timeStepWidth_ = this->specificSettings_.getOptionDouble("timestepWidth", 0.01, PythonUtility::Positive);
  LOG(DEBUG) << "precice initialization done, dt: " << maximumPreciceTimestepSize_ << "," << timeStepWidth_;

  initialized_ = true;

#else
  LOG(FATAL) << "Not compiled with preCICE!";
#endif
}

template<typename NestedSolver>
void MuscleContraction<NestedSolver>::
run()
{
#ifdef HAVE_PRECICE

  // initialize everything
  initialize();

  if (preciceSolverInterface_->isActionRequired(precice::constants::actionWriteInitialData()))
  {
    LOG(DEBUG) << "isActionRequired(actionWriteInitialData) is true";

    // writeData
    preciceWriteData();

    preciceSolverInterface_->markActionFulfilled(precice::constants::actionWriteInitialData());
    // initialize data in precice
    LOG(DEBUG) << "precice::initializeData";
    preciceSolverInterface_->initializeData();

  }
  else
  {
    LOG(DEBUG) << "isActionRequired(actionWriteInitialData) is false";
  }

  //std::shared_ptr<::Data::SlotConnectorData<::FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, ::BasisFunction::LagrangeOfOrder<2> >, 1, 1> > connectorData
  //  = getSlotConnectorData();

  // perform the computation of this solver

  // main simulation loop of adapter
  for (int timeStepNo = 0; preciceSolverInterface_->isCouplingOngoing(); timeStepNo++)
  {

    // read scalar gamma value from precice
    preciceReadData();

    // compute the time step width such that it fits in the remaining time in the current time window
    //int timeStepWidth = std::min(maximumPreciceTimestepSize_, timeStepWidth_);

    // hard-code 1 time step for the static problem
    double timeStepWidth = maximumPreciceTimestepSize_;

    // call the nested solver
    nestedSolver_.run();

    // write data to precice
    // data to send:
    // - geometry
    // - lambda (todo)
    // - lambdaDot (todo)
    preciceWriteData();

    LOG(DEBUG) << "precice::advance(" << timeStepWidth << "), maximumPreciceTimestepSize_: " << maximumPreciceTimestepSize_;

    // advance timestepping in precice
    maximumPreciceTimestepSize_ = preciceSolverInterface_->advance(timeStepWidth);
  }
  preciceSolverInterface_->finalize();

#endif
}

#ifdef HAVE_PRECICE

template<typename NestedSolver>
void MuscleContraction<NestedSolver>::
preciceReadData()
{
  if (!preciceSolverInterface_->isReadDataAvailable())
    return;

  LOG(DEBUG) << "read data from precice";

  std::shared_ptr<::FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, ::BasisFunction::LagrangeOfOrder<2>>> functionSpace
    = nestedSolver_.data().functionSpace();

  int nDofsLocalWithoutGhosts = functionSpace->nDofsLocalWithoutGhosts();

  // read scalar gamma value from precice
  std::vector<double> gammaValues(nDofsLocalWithoutGhosts);
  preciceSolverInterface_->readBlockScalarData(preciceDataIdGamma_, nDofsLocalWithoutGhosts, preciceVertexIds_.data(), gammaValues.data());

  // set gamma values
  getSlotConnectorData()->variable1[connectorSlotIdGamma_].setValuesWithoutGhosts(gammaValues);

  LOG(DEBUG) << "read data from precice complete, gamma values: " << gammaValues;
}

template<typename NestedSolver>
void MuscleContraction<NestedSolver>::
preciceWriteData()
{
  if (!preciceSolverInterface_->isWriteDataRequired(timeStepWidth_))
    return;

  std::shared_ptr<::FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, ::BasisFunction::LagrangeOfOrder<2>>> functionSpace
    = nestedSolver_.data().functionSpace();

  int nDofsLocalWithoutGhosts = functionSpace->nDofsLocalWithoutGhosts();

  // write data to precice

  // data to send:
  // - geometry
  // - lambda (todo)
  // - lambdaDot (todo)

  // convert geometry values to precice data layout
  std::vector<Vec3> geometryValues;
  functionSpace->geometryField().getValuesWithoutGhosts(geometryValues);

  std::vector<double> geometryValuesPrecice(3*nDofsLocalWithoutGhosts);

  for (dof_no_t dofNoLocal = 0; dofNoLocal < nDofsLocalWithoutGhosts; dofNoLocal++)
  {
    for (int i = 0; i < 3; i++)
    {
      geometryValuesPrecice[3*dofNoLocal + i] = geometryValues[dofNoLocal][i];
    }
  }

  LOG(DEBUG) << "write geometry data to precice: " << geometryValuesPrecice;

  // write geometry values in precice
  preciceSolverInterface_->writeBlockVectorData(preciceDataIdGeometry_, nDofsLocalWithoutGhosts,
                                                preciceVertexIds_.data(), geometryValuesPrecice.data());

  LOG(DEBUG) << "write geometry data to precice complete";
}

#endif

template<typename NestedSolver>
void MuscleContraction<NestedSolver>::
reset()
{
  nestedSolver_.reset();

  initialized_ = false;
  // "uninitialize" everything
}

template<typename NestedSolver>
typename MuscleContraction<NestedSolver>::Data &MuscleContraction<NestedSolver>::
data()
{
  // get a reference to the data object
  return nestedSolver_.data();
}

//! get the data that will be transferred in the operator splitting to the other term of the splitting
//! the transfer is done by the slot_connector_data_transfer class
template<typename NestedSolver>
std::shared_ptr<typename MuscleContraction<NestedSolver>::SlotConnectorDataType> MuscleContraction<NestedSolver>::
getSlotConnectorData()
{
  //! This is relevant only, if this solver is part of a splitting or coupling scheme. Then this method returns the values/variables that will be
  // transferred to the other solvers. We can just reuse the values of the nestedSolver_.
  return nestedSolver_.getSlotConnectorData();
}

}  // namespace
