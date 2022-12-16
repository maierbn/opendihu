#include "control/precice/surface_coupling/01_initialize.h"

#include <sstream>

namespace Control
{

template<typename NestedSolver>
PreciceAdapterInitialize<NestedSolver>::
PreciceAdapterInitialize(DihuContext context) :
  context_(context["PreciceAdapter"]),
  nestedSolver_(this->context_), maximumPreciceTimestepSize_(0), timeStepOutputInterval_(1), initialized_(false)
{
  // get python settings object from context
  this->specificSettings_ = this->context_.getPythonConfig();
}

template<typename NestedSolver>
void PreciceAdapterInitialize<NestedSolver>::
initialize()
{
#ifdef HAVE_PRECICE

  LOG(DEBUG) << "initialize precice adapter, initialized_=" << initialized_;

  // make sure that we initialize only once, in the next call, initialized_ is true
  if (initialized_)
    return;

  // initialize() will be called before the simulation starts.

  // add this solver to the solvers diagram, which is an ASCII art representation that will be created at the end of the simulation.
  DihuContext::solverStructureVisualizer()->addSolver("Control::PreciceAdapter", true);   // hasInternalConnectionToFirstNestedSolver=true (the last argument) means slot connector data is shared with the first subsolver

  // indicate in solverStructureVisualizer that now a child solver will be initialized
  DihuContext::solverStructureVisualizer()->beginChild();

  // call initialize of the nested solver
  nestedSolver_.initialize();

  // initialize function space
  functionSpace_ = this->functionSpace(nestedSolver_);

  // indicate in solverStructureVisualizer that the child solver initialization is done
  DihuContext::solverStructureVisualizer()->endChild();

  // set the slotConnectorData for the solverStructureVisualizer to appear in the solver diagram
  //DihuContext::solverStructureVisualizer()->setSlotConnectorData(this->getSlotConnectorData());

  // check if the coupling is enabled
  couplingEnabled_ = this->specificSettings_.getOptionBool("couplingEnabled", true);

  timeStepWidth_ = this->specificSettings_.getOptionDouble("timestepWidth", 0.01, PythonUtility::Positive);

  // if not enabled, abort initialization
  if (!couplingEnabled_)
  {
    endTimeIfCouplingDisabled_ = this->specificSettings_.getOptionDouble("endTimeIfCouplingDisabled", 1, PythonUtility::Positive);

    LOG(WARNING) << "Coupling in PreciceAdapterVolumeCoupling is disabled (option \"couplingEnabled\": False), "
      << "using end time \"endTimeIfCouplingDisabled\": " << endTimeIfCouplingDisabled_ << ".";

    initialized_ = true;
    return;
  }

  // initialize precice
  preciceParticipantName_ = this->specificSettings_.getOptionString("preciceParticipantName", "MuscleSolver");
  const std::string configFileName = this->specificSettings_.getOptionString("preciceConfigFilename", "../precice-config.xml");
  outputOnlyConvergedTimeSteps_ = this->specificSettings_.getOptionBool("outputOnlyConvergedTimeSteps", true);

  int rankNo = functionSpace_->meshPartition()->rankSubset()->ownRankNo();
  int nRanks = functionSpace_->meshPartition()->rankSubset()->size();

  // initialize interface to precice for the bottom surface mesh
  preciceSolverInterface_ = std::make_shared<precice::SolverInterface>(preciceParticipantName_, configFileName, rankNo, nRanks);

  // parse the options in "preciceMeshes" and initialize all meshes in precice, store in variable preciceMeshes_
  initializePreciceMeshes();

  // parse the options in "preciceData" and initialize all variables in precice, store in variable preciceData_
  initializePreciceData();

  // initialize Dirichlet boundary conditions at all dofs that will get some prescribed values during coupling
  initializeDirichletBoundaryConditions();

  // parse scalingFactor from settings
  scalingFactor_ = this->specificSettings_.getOptionDouble("scalingFactor", 1);

  // determine maximum timestep size
  maximumPreciceTimestepSize_ = std::max(maximumPreciceTimestepSize_, preciceSolverInterface_->initialize());

  LOG(DEBUG) << "precice initialization done, dt: " << maximumPreciceTimestepSize_ << "," << timeStepWidth_;

  initialized_ = true;

#else
  LOG(FATAL) << "Failed to initialize PreciceAdapter (surface coupling) because opendihu is not compiled with preCICE.";
#endif
}


#ifdef HAVE_PRECICE
template<typename NestedSolver>
void PreciceAdapterInitialize<NestedSolver>::
initializePreciceMeshes()
{
  const int nNodesX = functionSpace_->nNodesLocalWithoutGhosts(0);
  const int nNodesY = functionSpace_->nNodesLocalWithoutGhosts(1);
  const int nNodesZ = functionSpace_->nNodesLocalWithoutGhosts(2);

  std::vector<Vec3> geometryValues;
  functionSpace_->geometryField().getValuesWithoutGhosts(geometryValues);

  // parse settings of meshes
  std::string settingsKey("preciceMeshes");
  PyObject *listPy = this->specificSettings_.getOptionPyObject(settingsKey);
  std::vector<PyObject *> list = PythonUtility::convertFromPython<std::vector<PyObject *>>::get(listPy);
  PythonConfig preciceMeshConfig(this->specificSettings_, settingsKey);

  // loop over items of the list under "preciceMeshes"
  for (int i = 0; i < list.size(); i++)
  {
    PythonConfig currentMeshConfig(preciceMeshConfig, i);

    std::shared_ptr<PreciceMesh> preciceMesh = std::make_shared<PreciceMesh>();

    // parse name of mesh
    preciceMesh->preciceMeshName = currentMeshConfig.getOptionString("preciceMeshName", "");

    // get precice id from precice config
    preciceMesh->preciceMeshId = preciceSolverInterface_->getMeshID(preciceMesh->preciceMeshName);

    // parse face
    std::string face = currentMeshConfig.getOptionString("face", "2-");
    if (face == "2-")
    {
      preciceMesh->face = PreciceMesh::face2Minus;
    }
    else if (face == "2+")
    {
      preciceMesh->face = PreciceMesh::face2Plus;
    }
    else
    {
      LOG(FATAL) << currentMeshConfig << "[\"face\"] is \"" << face << "\", valid values are: \"2-\", \"2+\".";
    }

    // check if there are any local nodes of the surface on the local partition
    bool localDomainHasPartOfSurface = true;
    if (preciceMesh->face == PreciceMesh::face2Minus && functionSpace_->meshPartition()->ownRankPartitioningIndex(2) > 0)
    {
      localDomainHasPartOfSurface = false;
    }
    else if (preciceMesh->face == PreciceMesh::face2Plus && functionSpace_->meshPartition()->ownRankPartitioningIndex(2) < functionSpace_->meshPartition()->nRanks(2)-1)
    {
      localDomainHasPartOfSurface = false;
    }

    if (localDomainHasPartOfSurface)
    {
      // store number of nodes
      preciceMesh->nNodesLocal = nNodesX * nNodesY;

      // collect node positions for all surface nodes of the coupling surface
      std::vector<double> geometryValuesSurface(3*preciceMesh->nNodesLocal);

      // resize buffer for the local dof nos in the 3D mesh of the surface mesh
      preciceMesh->dofNosLocal.resize(preciceMesh->nNodesLocal);

      int nodeIndexZ = 0;
      if (face == "2+")
      {
        nodeIndexZ = (nNodesZ-1);
      }

      // loop over nodes
      for (int nodeIndexY = 0; nodeIndexY < nNodesY; nodeIndexY++)
      {
        for (int nodeIndexX = 0; nodeIndexX < nNodesX; nodeIndexX++)
        {
          node_no_t surfaceDofNo = nodeIndexY * nNodesX + nodeIndexX;
          dof_no_t dofNoLocal =
            nodeIndexZ * nNodesX * nNodesY
            + nodeIndexY * nNodesX
            + nodeIndexX;

          preciceMesh->dofNosLocal[surfaceDofNo] = dofNoLocal;

          for (int i = 0; i < 3; i++)
          {
            geometryValuesSurface[3*surfaceDofNo + i] = geometryValues[dofNoLocal][i];
          }
        }
      }

      // resize buffer for vertex ids
      preciceMesh->preciceVertexIds.resize(preciceMesh->nNodesLocal);

      // give the node positions to precice and get the vertex ids
      // void precice::SolverInterface::setMeshVertices(int meshID, int size, const double *positions, int *ids)
      preciceSolverInterface_->setMeshVertices(preciceMesh->preciceMeshId, preciceMesh->nNodesLocal, geometryValuesSurface.data(), preciceMesh->preciceVertexIds.data());
    }
    else
    {
      // there are no local nodes of this surface
      preciceMesh->nNodesLocal = 0;
    }

    // store the precice mesh to the vector of meshes
    preciceMeshes_.push_back(preciceMesh);
  }
}

template<typename NestedSolver>
void PreciceAdapterInitialize<NestedSolver>::
initializePreciceData()
{
  // parse settings for coupling participants / tendons
  // loop over items of the key "preciceData"
  std::string settingsKey("preciceData");
  PyObject *listPy = this->specificSettings_.getOptionPyObject(settingsKey);
  std::vector<PyObject *> list = PythonUtility::convertFromPython<std::vector<PyObject *>>::get(listPy);
  PythonConfig preciceDataConfig(this->specificSettings_, settingsKey);

  // loop over items of the list under "preciceData"
  for (int i = 0; i < list.size(); i++)
  {
    PythonConfig currentPreciceData(preciceDataConfig, i);

    PreciceData preciceData;

    // parse the mesh
    std::string preciceMeshName = currentPreciceData.getOptionString("preciceMeshName", "");
    int preciceMeshId = preciceSolverInterface_->getMeshID(preciceMeshName);

    // find the mesh in the already parsed precice meshes
    typename std::vector<std::shared_ptr<PreciceMesh>>::iterator iter
      = std::find_if(preciceMeshes_.begin(), preciceMeshes_.end(), [&preciceMeshId](std::shared_ptr<PreciceMesh> preciceMesh)
    {
      return preciceMesh->preciceMeshId == preciceMeshId;
    });

    if (iter == preciceMeshes_.end())
    {
      std::stringstream s;
      for (std::shared_ptr<PreciceMesh> preciceMesh : preciceMeshes_)
        s << " \"" << preciceMesh->preciceMeshName << "\"";
      LOG(FATAL) << currentPreciceData << "[\"preciceMeshName\"] = \"" << preciceMeshName << "\" could not be found, available precice meshes: " << s.str();
    }
    preciceData.preciceMesh = *iter;

    // parse mode and variable names
    std::string mode = currentPreciceData.getOptionString("mode", "");
    if (mode == "read-displacements-velocities")
    {
      preciceData.ioType = PreciceData::ioRead;
      preciceData.boundaryConditionType = PreciceData::bcTypeDirichlet;

      // get precice names of the variables
      preciceData.displacementsName = currentPreciceData.getOptionString("displacementsName", "Displacement");
      // preciceData.velocitiesName = currentPreciceData.getOptionString("velocitiesName", "Velocity");

      // get precice data ids
      preciceData.preciceDataIdDisplacements = preciceSolverInterface_->getDataID(
        preciceData.displacementsName, preciceData.preciceMesh->preciceMeshId);

      // preciceData.preciceDataIdVelocities = preciceSolverInterface_->getDataID(
      //   preciceData.velocitiesName, preciceData.preciceMesh->preciceMeshId);
    }
    else if (mode == "read-traction")
    {
      preciceData.ioType = PreciceData::ioRead;
      preciceData.boundaryConditionType = PreciceData::bcTypeNeumann;

      // get precice names of the variables
      preciceData.tractionName = currentPreciceData.getOptionString("tractionName", "Traction");

      // get precice data ids
      preciceData.preciceDataIdTraction = preciceSolverInterface_->getDataID(
        preciceData.tractionName, preciceData.preciceMesh->preciceMeshId);
    }
    else if (mode == "write-displacements-velocities")
    {
      preciceData.ioType = PreciceData::ioWrite;

      // get precice names of the variables
      preciceData.displacementsName = currentPreciceData.getOptionString("displacementsName", "Displacement");
      // preciceData.velocitiesName = currentPreciceData.getOptionString("velocitiesName", "Velocity");

      // get precice data ids
      preciceData.preciceDataIdDisplacements = preciceSolverInterface_->getDataID(
        preciceData.displacementsName, preciceData.preciceMesh->preciceMeshId);

      // preciceData.preciceDataIdVelocities = preciceSolverInterface_->getDataID(
      //   preciceData.velocitiesName, preciceData.preciceMesh->preciceMeshId);
    }
    else if (mode == "write-traction")
    {
      preciceData.ioType = PreciceData::ioWrite;

      // get precice names of the variables
      preciceData.tractionName = currentPreciceData.getOptionString("tractionName", "Traction");

      // get precice data ids
      preciceData.preciceDataIdTraction = preciceSolverInterface_->getDataID(
        preciceData.tractionName, preciceData.preciceMesh->preciceMeshId);
    }
    else
    {
      LOG(FATAL) << currentPreciceData << "[\"mode\"] is \"" << mode << "\", "
        << "possible values are: \"read-displacements-velocities\", \"read-traction\", \"write-displacements-velocities\", \"write-traction\".";
    }

    // store preciceData to vector
    preciceData_.push_back(preciceData);
  }
}

template<typename NestedSolver>
void PreciceAdapterInitialize<NestedSolver>::
initializeDirichletBoundaryConditions()
{
  using ElementWithNodesType = typename SpatialDiscretization::DirichletBoundaryConditionsBase<typename PreciceAdapterNestedSolver<NestedSolver>::FunctionSpace,6>::ElementWithNodes;

  // initialize dirichlet boundary conditions, set all Dirichlet boundary condition values that will be needed later to vector (0,0,0)
  std::vector<ElementWithNodesType> dirichletBoundaryConditionElements;

  int nElementsX = functionSpace_->meshPartition()->nElementsLocal(0);
  int nElementsY = functionSpace_->meshPartition()->nElementsLocal(1);
  int nElementsZ = functionSpace_->meshPartition()->nElementsLocal(2);

  std::set<int> elementIndicesZ;

  // loop over precice field variables to be transferred, collect surface meshes at bottom or top
  for (PreciceData &preciceData : preciceData_)
  {
    if (preciceData.ioType == PreciceData::ioRead && preciceData.boundaryConditionType == PreciceData::bcTypeDirichlet)
    {
      if (preciceData.preciceMesh->face == PreciceMesh::face2Minus)
      {
        elementIndicesZ.insert(0);
      }
      else if (preciceData.preciceMesh->face == PreciceMesh::face2Plus)
      {
        elementIndicesZ.insert(nElementsZ-1);
      }
    }
  }

  // for either or both of top and bottom surface coupling mesh
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
            elementWithNodes.elementalDofIndex.insert(std::pair<int,VecD<6>>(elementalDofIndex, VecD<6>{0,0,0,0,0,0}));
          }
        }
        dirichletBoundaryConditionElements.push_back(elementWithNodes);
      }
    }
  }

  // add dirichlet bc values for all nodes that will get a prescribed value during coupling
  this->addDirichletBoundaryConditions(nestedSolver_, dirichletBoundaryConditionElements);
}
#endif

}  // namespace
