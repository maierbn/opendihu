#include "control/precice/contraction_dirichlet_boundary_conditions.h"

#include <sstream>

namespace PreciceAdapter
{

template<typename NestedSolver>
ContractionDirichletBoundaryConditions<NestedSolver>::
ContractionDirichletBoundaryConditions(DihuContext context) :
  Runnable(),
  context_(context["PreciceContractionDirichletBoundaryConditions"]),
  nestedSolver_(this->context_), maximumPreciceTimestepSize_(0), timeStepOutputInterval_(1),
  haveCouplingSurfaceBottom_(false), haveCouplingSurfaceTop_(false), initialized_(false)
{
  // get python settings object from context
  this->specificSettings_ = this->context_.getPythonConfig();
}

template<typename NestedSolver>
void ContractionDirichletBoundaryConditions<NestedSolver>::
initialize()
{
#ifdef HAVE_PRECICE

  LOG(DEBUG) << "initialize precice adapter for muscle, initialized_=" << initialized_;

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

  // initialize function space
  functionSpace_ = nestedSolver_.timeStepping2().data().functionSpace();

  // indicate in solverStructureVisualizer that the child solver initialization is done
  DihuContext::solverStructureVisualizer()->endChild();

  // set the outputConnectorData for the solverStructureVisualizer to appear in the solver diagram
  DihuContext::solverStructureVisualizer()->setOutputConnectorData(getOutputConnectorData());

  // initialize precice
  const std::string solverName = "MuscleSolver";
  const std::string configFileName = this->specificSettings_.getOptionString("preciceConfigFilename", "../precice-config.xml");

  int rankNo = functionSpace_->meshPartition()->rankSubset()->ownRankNo();
  int nRanks = functionSpace_->meshPartition()->rankSubset()->size();

  // parse settings for coupling participants / tendons
  // loop over items of the key "CouplingParticipants"
  std::string settingsKey("couplingParticipants");
  PyObject *listPy = this->specificSettings_.getOptionPyObject(settingsKey);
  std::vector<PyObject *> list = PythonUtility::convertFromPython<std::vector<PyObject *>>::get(listPy);
  PythonConfig couplingParticipantConfig(this->specificSettings_, settingsKey);

  // loop over items of the list under "couplingParticipants"
  for (int i = 0; i < list.size(); i++)
  {
    PythonConfig currentCouplingParticipantConfig(couplingParticipantConfig, i);

    CouplingParticipant couplingParticipant;

    // initialize if this participant (tendon) is attached at bottom or top of muscle mesh
    couplingParticipant.isCouplingSurfaceBottom = currentCouplingParticipantConfig.getOptionBool("isCouplingSurfaceBottom", true);

    // initialize interface to precice for the bottom surface mesh
    couplingParticipant.preciceSolverInterface = std::make_shared<precice::SolverInterface>(solverName, configFileName, rankNo, nRanks);

    couplingParticipants_.push_back(couplingParticipant);

    if (couplingParticipant.isCouplingSurfaceBottom)
    {
      haveCouplingSurfaceBottom_ = true;

      // get the mesh id of the bottom mesh
      std::string meshName = "MuscleMeshBottom";
      preciceMeshIdBottom_ = couplingParticipant.preciceSolverInterface->getMeshID(meshName);
    }
    else
    {
      haveCouplingSurfaceTop_ = true;

      // get the mesh id of the top mesh
      std::string meshName = "MuscleMeshTop";
      preciceMeshIdTop_ = couplingParticipant.preciceSolverInterface->getMeshID(meshName);
    }
  }

  initializeDirichletBoundaryConditions();

  std::vector<Vec3> geometryValues;
  functionSpace_->geometryField().getValuesWithoutGhosts(geometryValues);

  // store the node positions to precice

  // initialize coupling surfaces
  // get nodes at coupling surfaces
  const int nNodesX = functionSpace_->nNodesLocalWithoutGhosts(0);
  const int nNodesY = functionSpace_->nNodesLocalWithoutGhosts(1);
  const int nNodesZ = functionSpace_->nNodesLocalWithoutGhosts(2);
  nNodesSurfaceLocal_ = nNodesX * nNodesY;

  std::vector<double> geometryValuesBottomSurfacePrecice(3*nNodesSurfaceLocal_);
  std::vector<double> geometryValuesTopSurfacePrecice(3*nNodesSurfaceLocal_);

  // loop over nodes
  for (int nodeIndexY = 0; nodeIndexY < nNodesY; nodeIndexY++)
  {
    for (int nodeIndexX = 0; nodeIndexX < nNodesX; nodeIndexX++)
    {
      node_no_t nodeNoLocal = nodeIndexY * nNodesX + nodeIndexX;
      dof_no_t bottomDofNoLocal = nodeIndexY * nNodesX + nodeIndexX;
      dof_no_t topDofNoLocal =
        (nNodesZ-1) * nNodesX * nNodesY
        + nodeIndexY * nNodesX
        + nodeIndexX;

      for (int i = 0; i < 3; i++)
      {
        geometryValuesBottomSurfacePrecice[3*nodeNoLocal + i] = geometryValues[bottomDofNoLocal][i];
        geometryValuesTopSurfacePrecice[3*nodeNoLocal + i] = geometryValues[topDofNoLocal][i];
      }
    }
  }

  LOG(DEBUG) << "setMeshVertices to precice for surface meshes, " << 3*nNodesSurfaceLocal_ << " values each";

  if (haveCouplingSurfaceBottom_)
  {
    preciceVertexIdsBottom_.resize(nNodesSurfaceLocal_);

    // initialize all precice solver interfaces that use the bottom interface

    // loop over coupling participants
    for (CouplingParticipant &couplingParticipant : couplingParticipants_)
    {
      if (couplingParticipant.isCouplingSurfaceBottom)
      {
        std::shared_ptr<precice::SolverInterface> preciceSolverInterface = couplingParticipant.preciceSolverInterface;

        //preciceSolverInterface->setMeshVertices(preciceMeshId_, int size, double* positions, int* ids);
        preciceSolverInterface->setMeshVertices(preciceMeshIdBottom_, nNodesSurfaceLocal_, geometryValuesBottomSurfacePrecice.data(), preciceVertexIdsBottom_.data());

        LOG(DEBUG) << "precice defined vertexIds: " << preciceVertexIdsBottom_;

        // initialize data ids
        preciceDataIdDisplacements_ = preciceSolverInterface->getDataID("Displacements", preciceMeshIdBottom_);
        preciceDataIdVelocity_      = preciceSolverInterface->getDataID("Velocity",      preciceMeshIdBottom_);
        preciceDataIdTraction_      = preciceSolverInterface->getDataID("Traction",      preciceMeshIdBottom_);

        LOG(DEBUG) << "data id displacements: " << preciceDataIdDisplacements_;
        LOG(DEBUG) << "data id traction: " << preciceDataIdTraction_;

        maximumPreciceTimestepSize_ = std::max(maximumPreciceTimestepSize_, preciceSolverInterface->initialize());
      }
    }
  }

  if (haveCouplingSurfaceTop_)
  {
    preciceVertexIdsTop_.resize(nNodesSurfaceLocal_);

    // loop over coupling participants
    for (CouplingParticipant &couplingParticipant : couplingParticipants_)
    {
      if (!couplingParticipant.isCouplingSurfaceBottom)
      {
        std::shared_ptr<precice::SolverInterface> preciceSolverInterface = couplingParticipant.preciceSolverInterface;

        //preciceSolverInterface->setMeshVertices(preciceMeshId_, int size, double* positions, int* ids);
        preciceSolverInterface->setMeshVertices(preciceMeshIdTop_, nNodesSurfaceLocal_, geometryValuesTopSurfacePrecice.data(), preciceVertexIdsTop_.data());

        LOG(DEBUG) << "precice defined vertexIds: " << preciceVertexIdsTop_;

        // initialize data ids
        preciceDataIdDisplacements_ = preciceSolverInterface->getDataID("Displacements", preciceMeshIdTop_);
        preciceDataIdVelocity_      = preciceSolverInterface->getDataID("Velocity",      preciceMeshIdTop_);
        preciceDataIdTraction_      = preciceSolverInterface->getDataID("Traction",      preciceMeshIdTop_);

        LOG(DEBUG) << "data id displacements: " << preciceDataIdDisplacements_;
        LOG(DEBUG) << "data id traction: " << preciceDataIdTraction_;

        maximumPreciceTimestepSize_ = std::max(maximumPreciceTimestepSize_, preciceSolverInterface->initialize());
      }
    }
  }

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
  int nElementsZ = functionSpace_->meshPartition()->nElementsLocal(2);

  std::vector<int> elementIndicesZ;

  if (haveCouplingSurfaceBottom_)
  {
    elementIndicesZ.push_back(0);
  }

  if (haveCouplingSurfaceTop_)
  {
    elementIndicesZ.push_back(nElementsZ-1);
  }

  // for either or both of surface and bottom coupling mesh
  for (int elementIndexZ : elementIndicesZ)
  {
    int indexZ = 0;
    if (elementIndexZ > 0)
      indexZ = 2;

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
            int elementalDofIndex = indexZ * 9 + indexY * 3 + indexX;
            elementWithNodes.elementalDofIndex.push_back(std::pair<int,VecD<6>>(elementalDofIndex, VecD<6>{0,0,0,0,0,0}));
          }
        }
        dirichletBoundaryConditionElements.push_back(elementWithNodes);
      }
    }
  }

  // add the dirichlet bc values
  bool overwriteBcOnSameDof = true;
  nestedSolver_.timeStepping2().dynamicHyperelasticitySolver()->addDirichletBoundaryConditions(dirichletBoundaryConditionElements, overwriteBcOnSameDof);
}

template<typename NestedSolver>
void ContractionDirichletBoundaryConditions<NestedSolver>::
run()
{
#ifdef HAVE_PRECICE

  // initialize everything
  initialize();

  // perform initial data transfer for all participants, if required
  for (const CouplingParticipant &couplingParticipant : couplingParticipants_)
  {
    if (couplingParticipant.preciceSolverInterface->isActionRequired(precice::constants::actionWriteInitialData()))
    {
      // writeData for this participant
      preciceWriteData(couplingParticipant);

      couplingParticipant.preciceSolverInterface->markActionFulfilled(precice::constants::actionWriteInitialData());

      // initialize data in precice
      couplingParticipant.preciceSolverInterface->initializeData();
    }
  }

  // perform the computation of this solver
  double currentTime = 0;

  // main simulation loop of adapter
  for (int timeStepNo = 0; couplingParticipants_[0].preciceSolverInterface->isCouplingOngoing(); timeStepNo++)
  {
    if (timeStepNo % this->timeStepOutputInterval_ == 0 && (this->timeStepOutputInterval_ <= 10 || timeStepNo > 0))  // show first timestep only if timeStepOutputInterval is <= 10
    {
      LOG(INFO) << "Precice (Dirichlet) coupling, timestep " << timeStepNo << ", t=" << currentTime;
    }
    // read displacement values
    for (const CouplingParticipant &couplingParticipant : couplingParticipants_)
    {
      preciceReadData(couplingParticipant);
    }

    // compute the time step width such that it fits in the remaining time in the current time window
    double timeStepWidth = std::min(maximumPreciceTimestepSize_, timeStepWidth_);

    // set time span in nested solver
    nestedSolver_.setTimeSpan(currentTime, currentTime+timeStepWidth);

    // call the nested solver to proceed with the simulation for the assigned time span
    nestedSolver_.advanceTimeSpan();

    // write traction data to precice
    for (const CouplingParticipant &couplingParticipant : couplingParticipants_)
    {
      preciceWriteData(couplingParticipant);
    }

    // increase current simulation time
    currentTime += timeStepWidth;

    LOG(DEBUG) << "precice::advance(" << timeStepWidth << "), maximumPreciceTimestepSize_: " << maximumPreciceTimestepSize_;

    // advance timestepping in precice
    maximumPreciceTimestepSize_ = 0;
    for (const CouplingParticipant &couplingParticipant : couplingParticipants_)
    {
      double maximumPreciceTimestepSize = couplingParticipant.preciceSolverInterface->advance(timeStepWidth);
      maximumPreciceTimestepSize_ = std::max(maximumPreciceTimestepSize_, maximumPreciceTimestepSize);
    }
  }   // loop over time steps

  // finalize all precice interfaces
  for (const CouplingParticipant &couplingParticipant : couplingParticipants_)
  {
    couplingParticipant.preciceSolverInterface->finalize();
  }

#else
  LOG(FATAL) << "Not compiled with preCICE!";
#endif
}

#ifdef HAVE_PRECICE

template<typename NestedSolver>
void ContractionDirichletBoundaryConditions<NestedSolver>::
preciceReadData(const ContractionDirichletBoundaryConditions<NestedSolver>::CouplingParticipant &couplingParticipant)
{
  if (!couplingParticipant.preciceSolverInterface->isReadDataAvailable())
    return;

  LOG(DEBUG) << "read data from precice (displacement)";

  // bottom coupling surface

  // read displacement values from precice
  std::vector<double> displacementValues(nNodesSurfaceLocal_*3);
  std::vector<double> velocityValues(nNodesSurfaceLocal_*3);

  if (couplingParticipant.isCouplingSurfaceBottom)
  {
    couplingParticipant.preciceSolverInterface->readBlockVectorData(preciceDataIdDisplacements_, nNodesSurfaceLocal_,
                                                                    preciceVertexIdsBottom_.data(), displacementValues.data());

    couplingParticipant.preciceSolverInterface->readBlockVectorData(preciceDataIdVelocity_, nNodesSurfaceLocal_,
                                                                    preciceVertexIdsBottom_.data(), velocityValues.data());
  }
  else
  {
    couplingParticipant.preciceSolverInterface->readBlockVectorData(preciceDataIdDisplacements_, nNodesSurfaceLocal_,
                                                                    preciceVertexIdsTop_.data(), displacementValues.data());

    couplingParticipant.preciceSolverInterface->readBlockVectorData(preciceDataIdVelocity_, nNodesSurfaceLocal_,
                                                                    preciceVertexIdsTop_.data(), velocityValues.data());
  }
  // loop over nodes to set the received values
  std::vector<std::pair<global_no_t,std::array<double,6>>> newDirichletBCValues;
  newDirichletBCValues.reserve(nNodesSurfaceLocal_);

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

  // The nested solver is Control::Coupling<   (electrophysiology), MuscleContractionSolver<> >
  // Set the bc in the dynamic solver of the MuscleContractionSolver
  nestedSolver_.timeStepping2().dynamicHyperelasticitySolver()->updateDirichletBoundaryConditions(newDirichletBCValues);
}

template<typename NestedSolver>
void ContractionDirichletBoundaryConditions<NestedSolver>::
preciceWriteData(const ContractionDirichletBoundaryConditions<NestedSolver>::CouplingParticipant &couplingParticipant)
{
  if (!couplingParticipant.preciceSolverInterface->isWriteDataRequired(timeStepWidth_))
    return;

  // write traction data to precice

  // convert geometry values to precice data layout
  std::vector<Vec3> tractionValues;
  nestedSolver_.timeStepping2().data().materialTraction()->getValuesWithoutGhosts(tractionValues);

  std::vector<double> tractionValuesPrecice(3*nNodesSurfaceLocal_);

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

  if (couplingParticipant.isCouplingSurfaceBottom)
  {
    // write geometry values in precice
    couplingParticipant.preciceSolverInterface->writeBlockVectorData(preciceDataIdTraction_, nNodesSurfaceLocal_,
                                                                     preciceVertexIdsBottom_.data(), tractionValuesPrecice.data());
  }
  else
  {
    // write geometry values in precice
    couplingParticipant.preciceSolverInterface->writeBlockVectorData(preciceDataIdTraction_, nNodesSurfaceLocal_,
                                                                     preciceVertexIdsTop_.data(), tractionValuesPrecice.data());
  }

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
