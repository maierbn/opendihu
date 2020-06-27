#include "control/precice/contraction_dirichlet_boundary_conditions.h"

#include <sstream>

namespace PreciceAdapter
{

template<typename NestedSolver>
ContractionDirichletBoundaryConditions<NestedSolver>::
ContractionDirichletBoundaryConditions(DihuContext context) :
  Runnable(),
  context_(context["PreciceContractionDirichletBoundaryConditions"]), nestedSolver_(this->context_), initialized_(false)
{
  // get python settings object from context
  this->specificSettings_ = this->context_.getPythonConfig();
}

template<typename NestedSolver>
void ContractionDirichletBoundaryConditions<NestedSolver>::
initialize()
{
#ifdef HAVE_PRECICE

  // make sure that we initialize only once, in the next call, initialized_ is true
  if (initialized_)
    return;

  // initialize() will be called before the simulation starts.

  // add this solver to the solvers diagram, which is an ASCII art representation that will be created at the end of the simulation.
  DihuContext::solverStructureVisualizer()->addSolver("PreciceAdapter::ContractionDirichletBoundaryConditions", true);   // hasInternalConnectionToFirstNestedSolver=true (the last argument) means output connector data is shared with the first subsolver

  // indicate in solverStructureVisualizer that now a child solver will be initialized
  DihuContext::solverStructureVisualizer()->beginChild();

  // call initialize of the nested solver
  nestedSolver_.initialize();

  initializeDirichletBoundaryConditions();

  // indicate in solverStructureVisualizer that the child solver initialization is done
  DihuContext::solverStructureVisualizer()->endChild();

  // set the outputConnectorData for the solverStructureVisualizer to appear in the solver diagram
  DihuContext::solverStructureVisualizer()->setOutputConnectorData(getOutputConnectorData());

  // initialize precice
  const std::string solverName = "MuscleSolver";
  const std::string configFileName = this->specificSettings_.getOptionString("preciceConfigFilename", "../precice-config.xml");

  // initialize function space
  functionSpace_ = nestedSolver_.timeStepping2().data().functionSpace();

  int rankNo = functionSpace_->meshPartition()->rankSubset()->ownRankNo();
  int nRanks = functionSpace_->meshPartition()->rankSubset()->size();

  // initialize interface to precice for the bottom surface mesh
  preciceSolverInterfaceBottom_ = std::make_unique<precice::SolverInterface>(solverName, configFileName, rankNo, nRanks);

  std::vector<Vec3> geometryValues;
  functionSpace_->geometryField().getValuesWithoutGhosts(geometryValues);

  // store the node positions to precice
  std::string meshName = "MuscleMeshBottom";
  preciceMeshIdBottom_ = preciceSolverInterfaceBottom_->getMeshID(meshName);

  // bottom coupling surface
  // get nodes at coupling surface
  const int nNodesX = functionSpace_->nNodesLocalWithoutGhosts(0);
  const int nNodesY = functionSpace_->nNodesLocalWithoutGhosts(1);
  nNodesBottomSurfaceLocal_ = nNodesX * nNodesY;

  std::vector<double> geometryValuesSurfacePrecice(3*nNodesBottomSurfaceLocal_);
  preciceVertexIdsBottom_.resize(nNodesBottomSurfaceLocal_);

  // loop over nodes
  const int nodeIndexZ = 0;
  for (int nodeIndexY = 0; nodeIndexY < nNodesY; nodeIndexY++)
  {
    for (int nodeIndexX = 0; nodeIndexX < nNodesX; nodeIndexX++)
    {
      node_no_t nodeNoLocal =
        nodeIndexZ * nNodesX * nNodesY
        + nodeIndexY * nNodesX
        + nodeIndexX;
      dof_no_t dofNoLocal = nodeNoLocal;

      for (int i = 0; i < 3; i++)
      {
        geometryValuesSurfacePrecice[3*nodeNoLocal + i] = geometryValues[dofNoLocal][i];
      }
    }
  }

  LOG(DEBUG) << "setMeshVertices to precice for bottom mesh, " << 3*nNodesBottomSurfaceLocal_ << " values";

  //preciceSolverInterfaceBottom_->setMeshVertices(preciceMeshId_, int size, double* positions, int* ids);
  preciceSolverInterfaceBottom_->setMeshVertices(preciceMeshIdBottom_, nNodesBottomSurfaceLocal_, geometryValuesSurfacePrecice.data(), preciceVertexIdsBottom_.data());

  LOG(DEBUG) << "precice defined vertexIds: " << preciceVertexIdsBottom_;

  // initialize data ids
  preciceDataIdDisplacements_ = preciceSolverInterfaceBottom_->getDataID("Displacements", preciceMeshIdBottom_);
  preciceDataIdTraction_      = preciceSolverInterfaceBottom_->getDataID("Traction",      preciceMeshIdBottom_);

  LOG(DEBUG) << "data id displacements: " << preciceDataIdDisplacements_;
  LOG(DEBUG) << "data id traction: " << preciceDataIdTraction_;

  maximumPreciceTimestepSize_ = preciceSolverInterfaceBottom_->initialize();

  timeStepWidth_ = this->specificSettings_.getOptionDouble("timestepWidth", 0.01, PythonUtility::Positive);
  LOG(DEBUG) << "precice initialization done, dt: " << maximumPreciceTimestepSize_ << "," << timeStepWidth_;


  initialized_ = true;

#else
  LOG(FATAL) << "Not compiled with preCICE!";
#endif
}

template<typename NestedSolver>
void ContractionDirichletBoundaryConditions<NestedSolver>::
initializeDirichletBoundaryConditions()
{
  using ElementWithNodesType = typename SpatialDiscretization::DirichletBoundaryConditionsBase<FunctionSpace,6>::ElementWithNodes;

  // initialize dirichlet boundary conditions, set all Dirichlet boundary condition values that will be needed later to vector (0,0,0)
  std::vector<ElementWithNodesType> dirichletBoundaryConditionElements;

  int nElementsX = functionSpace_->meshPartition()->nElementsLocal(0);
  int nElementsY = functionSpace_->meshPartition()->nElementsLocal(1);
  int elementIndexZ = 0;

  // loop over elements
  for (int elementIndexY = 0; elementIndexY < nElementsY; elementIndexY++)
  {
    for (int elementIndexX = 0; elementIndexX < nElementsX; elementIndexX++)
    {
      ElementWithNodesType elementWithNodes;
      elementWithNodes.elementNoLocal = elementIndexZ*nElementsX*nElementsY + elementIndexY * nElementsX + elementIndexX;

      for (int indexY = 0; indexY < 3; indexY++)
      {
        for (int indexX = 0; indexX < 3; indexX++)
        {
          int elementalDofIndex = indexY * 3 + indexX;
          elementWithNodes.elementalDofIndex.push_back(std::pair<int,VecD<6>>(elementalDofIndex, VecD<6>{0,0,0,0,0,0}));
        }
      }
      dirichletBoundaryConditionElements.push_back(elementWithNodes);
    }
  }

  // add the dirichlet bc values
  bool overwriteBcOnSameDof = true;
  nestedSolver_.timeStepping2().dynamicHyperelasticitySolver()->hyperelasticitySolver().dirichletBoundaryConditions()->addBoundaryConditions(dirichletBoundaryConditionElements, overwriteBcOnSameDof);
}

template<typename NestedSolver>
void ContractionDirichletBoundaryConditions<NestedSolver>::
run()
{
#ifdef HAVE_PRECICE

  // initialize everything
  initialize();

  if (preciceSolverInterfaceBottom_->isActionRequired(precice::constants::actionWriteInitialData()))
  {
    LOG(DEBUG) << "isActionRequired(actionWriteInitialData) is true";

    // writeData
    preciceWriteData();

    preciceSolverInterfaceBottom_->markActionFulfilled(precice::constants::actionWriteInitialData());
    // initialize data in precice
    LOG(DEBUG) << "precice::initializeData";
    preciceSolverInterfaceBottom_->initializeData();

  }
  else
  {
    LOG(DEBUG) << "isActionRequired(actionWriteInitialData) is false";
  }

  //std::shared_ptr<::Data::OutputConnectorData<::FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, ::BasisFunction::LagrangeOfOrder<2> >, 1, 1> > connectorData
  //  = getOutputConnectorData();

  // perform the computation of this solver

  // main simulation loop of adapter
  for (int timeStepNo = 0; preciceSolverInterfaceBottom_->isCouplingOngoing(); timeStepNo++)
  {

    // read displacement values
    preciceReadData();

    // compute the time step width such that it fits in the remaining time in the current time window
    //int timeStepWidth = std::min(maximumPreciceTimestepSize_, timeStepWidth_);

    // hard-code 1 time step for the static problem
    double timeStepWidth = maximumPreciceTimestepSize_;

    // call the nested solver
    nestedSolver_.run();

    // write traction data to precice
    preciceWriteData();

    LOG(DEBUG) << "precice::advance(" << timeStepWidth << "), maximumPreciceTimestepSize_: " << maximumPreciceTimestepSize_;

    // advance timestepping in precice
    maximumPreciceTimestepSize_ = preciceSolverInterfaceBottom_->advance(timeStepWidth);
  }
  preciceSolverInterfaceBottom_->finalize();

#endif
}

#ifdef HAVE_PRECICE

template<typename NestedSolver>
void ContractionDirichletBoundaryConditions<NestedSolver>::
preciceReadData()
{
  if (!preciceSolverInterfaceBottom_->isReadDataAvailable())
    return;

  LOG(DEBUG) << "read data from precice (displacement)";

  // bottom coupling surface

  // read displacement values from precice
  std::vector<double> displacementValues(nNodesBottomSurfaceLocal_*3);
  std::vector<double> velocityValues(nNodesBottomSurfaceLocal_*3);

  preciceSolverInterfaceBottom_->readBlockVectorData(preciceDataIdDisplacements_, nNodesBottomSurfaceLocal_,
                                                     preciceVertexIdsBottom_.data(), displacementValues.data());

  preciceSolverInterfaceBottom_->readBlockVectorData(preciceDataIdVelocity_, nNodesBottomSurfaceLocal_,
                                                     preciceVertexIdsBottom_.data(), velocityValues.data());
  // loop over nodes to set the received values
  std::vector<std::pair<global_no_t,std::array<double,6>>> newDirichletBCValues;
  newDirichletBCValues.reserve(nNodesBottomSurfaceLocal_);

  std::vector<double> geometryValuesSurfacePrecice(3*nNodesBottomSurfaceLocal_);
  preciceVertexIdsBottom_.resize(nNodesBottomSurfaceLocal_);

  // loop over nodes
  const int nNodesX = functionSpace_->nNodesLocalWithoutGhosts(0);
  const int nNodesY = functionSpace_->nNodesLocalWithoutGhosts(1);
  const int nodeIndexZ = 0;
  for (int nodeIndexY = 0; nodeIndexY < nNodesY; nodeIndexY++)
  {
    for (int nodeIndexX = 0; nodeIndexX < nNodesX; nodeIndexX++)
    {
      node_no_t nodeNoLocal =
        nodeIndexZ * nNodesX * nNodesY
        + nodeIndexY * nNodesX
        + nodeIndexX;

      dof_no_t dofNoLocal = nodeNoLocal;
      global_no_t dofNoGlobal = functionSpace_->meshPartition()->getDofNoGlobalPetsc(dofNoLocal);

      // assign received values to dirichlet bc vector of size 6
      std::array<double,6> newDirichletBCValue;

      for (int i = 0; i < 3; i++)
      {
        newDirichletBCValue[i] = displacementValues[3*dofNoLocal + i];
        newDirichletBCValue[3+i] = velocityValues[3*dofNoLocal + i];
      }

      newDirichletBCValues.push_back(std::pair<global_no_t,std::array<double,6>>(dofNoGlobal, newDirichletBCValue));
    }
  }

  LOG(DEBUG) << "read data from precice complete, displacement values: " << displacementValues << ", velocityValues: " << velocityValues;
  LOG(DEBUG) << "dirichlet bc to set: " << newDirichletBCValues;

  //! set new dirichlet boundary condition values
  nestedSolver_.timeStepping2().dynamicHyperelasticitySolver()->updateDirichletBoundaryConditions(newDirichletBCValues);
}

template<typename NestedSolver>
void ContractionDirichletBoundaryConditions<NestedSolver>::
preciceWriteData()
{
  if (!preciceSolverInterfaceBottom_->isWriteDataRequired(timeStepWidth_))
    return;

  // write traction data to precice

  // convert geometry values to precice data layout
  std::vector<Vec3> tractionValues;
  nestedSolver_.timeStepping2().data().materialTraction()->getValuesWithoutGhosts(tractionValues);

  std::vector<double> tractionValuesPrecice(3*nNodesBottomSurfaceLocal_);

  const int nNodesX = functionSpace_->nNodesLocalWithoutGhosts(0);
  const int nNodesY = functionSpace_->nNodesLocalWithoutGhosts(1);
  const int nodeIndexZ = 0;
  for (int nodeIndexY = 0; nodeIndexY < nNodesY; nodeIndexY++)
  {
    for (int nodeIndexX = 0; nodeIndexX < nNodesX; nodeIndexX++)
    {
      node_no_t nodeNoLocal =
        nodeIndexZ * nNodesX * nNodesY
        + nodeIndexY * nNodesX
        + nodeIndexX;

      dof_no_t dofNoLocal = nodeNoLocal;
      for (int i = 0; i < 3; i++)
      {
        tractionValuesPrecice[3*dofNoLocal + i] = tractionValues[dofNoLocal][i];
      }
    }
  }

  LOG(DEBUG) << "write traction data to precice: " << tractionValuesPrecice;

  // write geometry values in precice
  preciceSolverInterfaceBottom_->writeBlockVectorData(preciceDataIdTraction_, nNodesBottomSurfaceLocal_,
                                                      preciceVertexIdsBottom_.data(), tractionValuesPrecice.data());

  LOG(DEBUG) << "write traction data to precice complete";
}

#endif

template<typename NestedSolver>
void ContractionDirichletBoundaryConditions<NestedSolver>::
reset()
{
  nestedSolver_.reset();

  initialized_ = false;
  // "uninitialize" everything
}

template<typename NestedSolver>
typename ContractionDirichletBoundaryConditions<NestedSolver>::Data &ContractionDirichletBoundaryConditions<NestedSolver>::
data()
{
  // get a reference to the data object
  return nestedSolver_.data();
}

//! get the data that will be transferred in the operator splitting to the other term of the splitting
//! the transfer is done by the output_connector_data_transfer class
template<typename NestedSolver>
std::shared_ptr<typename ContractionDirichletBoundaryConditions<NestedSolver>::OutputConnectorDataType> ContractionDirichletBoundaryConditions<NestedSolver>::
getOutputConnectorData()
{
  //! This is relevant only, if this solver is part of a splitting or coupling scheme. Then this method returns the values/variables that will be
  // transferred to the other solvers. We can just reuse the values of the nestedSolver_.
  return nestedSolver_.getOutputConnectorData();
}

}  // namespace
