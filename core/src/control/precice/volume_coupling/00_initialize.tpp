#include "control/precice/volume_coupling/00_initialize.h"

#include <sstream>

#include "slot_connection/data_helper/slot_connector_data_helper.h"
#include "partition/mesh_partition/00_mesh_partition_base.h"
#include "mesh/mesh_manager/mesh_manager.h"

namespace Control
{

template<typename NestedSolver>
PreciceAdapterVolumeCouplingInitialize<NestedSolver>::
PreciceAdapterVolumeCouplingInitialize(DihuContext context) :
  context_(context["PreciceAdapterVolumeCoupling"]),
  nestedSolver_(this->context_), maximumPreciceTimestepSize_(0), timeStepOutputInterval_(1), initialized_(false)
{
  // get python settings object from context
  this->specificSettings_ = this->context_.getPythonConfig();
}

template<typename NestedSolver>
void PreciceAdapterVolumeCouplingInitialize<NestedSolver>::
initialize()
{
#ifdef HAVE_PRECICE

  LOG(DEBUG) << "initialize precice adapter (volume coupling), initialized_=" << initialized_;

  // make sure that we initialize only once, in the next call, initialized_ is true
  if (initialized_)
    return;

  // initialize() will be called before the simulation starts.

  // add this solver to the solvers diagram, which is an ASCII art representation that will be created at the end of the simulation.
  DihuContext::solverStructureVisualizer()->addSolver("Control::PreciceAdapterVolumeCoupling", true);   // hasInternalConnectionToFirstNestedSolver=true (the last argument) means slot connector data is shared with the first subsolver

  // indicate in solverStructureVisualizer that now a child solver will be initialized
  DihuContext::solverStructureVisualizer()->beginChild();

  // call initialize of the nested solver
  nestedSolver_.initialize();

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
  const std::string configFileName = this->specificSettings_.getOptionString("preciceConfigFilename", "../precice_config.xml");
  outputOnlyConvergedTimeSteps_ = this->specificSettings_.getOptionBool("outputOnlyConvergedTimeSteps", true);

  // get the union of all MPI ranks that occur for any fiber
  std::shared_ptr<Partition::RankSubset> rankSubset = this->context_.partitionManager()->rankSubsetForCollectiveOperations();
  int rankNo = rankSubset->ownRankNo();
  int nRanks = rankSubset->size();

  // debugging output of all ranks
  //std::cout << "rank for precice: " << rankNo << "/" << nRanks << " " << DihuContext::ownRankNoCommWorld() << "/" << DihuContext::nRanksCommWorld() << std::endl;

  // initialize interface to precice 
  preciceSolverInterface_ = std::make_shared<precice::SolverInterface>(preciceParticipantName_, configFileName, rankNo, nRanks);

  // parse the options in "preciceData" and initialize all variables in precice, store in variable preciceData_
  initializePreciceData();

  // parse scalingFactor from settings
  scalingFactor_ = this->specificSettings_.getOptionDouble("scalingFactor", 1);

  // determine maximum timestep size
  maximumPreciceTimestepSize_ = std::max(maximumPreciceTimestepSize_, preciceSolverInterface_->initialize());
  LOG(DEBUG) << "precice initialization done, dt: " << maximumPreciceTimestepSize_ << "," << timeStepWidth_;

  initialized_ = true;

#else
  LOG(FATAL) << "Failed to initialize PreciceAdapterVolumeCoupling because opendihu is not compiled with preCICE.";
#endif
}


//! get the nested solver
template<typename NestedSolver>
NestedSolver &PreciceAdapterVolumeCouplingInitialize<NestedSolver>::
nestedSolver()
{
  return nestedSolver_;
}

#ifdef HAVE_PRECICE

template<typename NestedSolver>
void PreciceAdapterVolumeCouplingInitialize<NestedSolver>::
initializePreciceData()
{
  // parse settings for coupling participants
  // loop over items of the key "preciceData"
  std::string settingsKey("preciceData");
  PyObject *listPy = this->specificSettings_.getOptionPyObject(settingsKey);
  std::vector<PyObject *> list = PythonUtility::convertFromPython<std::vector<PyObject *>>::get(listPy);
  PythonConfig preciceDataConfig(this->specificSettings_, settingsKey);

  // get the slot connector data from the nested solver
  int slotNo = 0;
  using SlotConnectorDataType = typename NestedSolver::SlotConnectorDataType;
  std::shared_ptr<SlotConnectorDataType> slotConnectorData = nestedSolver_.getSlotConnectorData();

  // get all present slot names
  std::vector<std::string> slotNames;
  SlotConnectorDataHelper<SlotConnectorDataType>::getSlotNames(slotConnectorData, slotNames);

  // loop over items of the list under "preciceData", the order of the items is also the order of the data connector slots
  for (int listIndex = 0; listIndex < list.size(); listIndex++)
  {
    // extract the list item as PythonConfig {"mode": "write", "variableName": ...}
    PythonConfig currentPreciceData(preciceDataConfig, listIndex);

    // create a new preciceData instance that will be added to preciceData_ at the end of this loop
    PreciceData preciceData;

    // parse options
    preciceData.isGeometryField = currentPreciceData.getOptionBool("isGeometryField", false);
    preciceData.opendihuMeshName = currentPreciceData.getOptionString("opendihuMeshName", "");

    // parse slot name
    preciceData.slotName = currentPreciceData.getOptionString("slotName", "");
    preciceData.slotNo = slotNo;

    LOG(DEBUG) << "isGeometryField: " << preciceData.isGeometryField << ", slotName: \"" << preciceData.slotName << "\".";

    if (preciceData.slotName == "")
    {
      preciceData.slotNo = slotNo;
      LOG(DEBUG) << "Using slotNo " << slotNo << " because slotName was not given.";
    }
    else
    {
      LOG(DEBUG) << "The slotName \"" << preciceData.slotName << "\" is given.";
      std::vector<std::string>::iterator slotNameIter = std::find(slotNames.begin(), slotNames.end(), preciceData.slotName);

      // if the "from" and "to" slot names match to slot names in the current slotConnectorData
      if (slotNameIter != slotNames.end())
      {
        LOG(DEBUG) << "Slot found";
        // determine slot nos
        preciceData.slotNo = std::distance(slotNames.begin(), slotNameIter);
        LOG(DEBUG) << "Slot \"" << preciceData.slotName << "\" is slot no " << preciceData.slotNo;
      }
      else
      {
        LOG(DEBUG) << "Slot not found";
        preciceData.slotNo = slotNo;
        LOG(WARNING) << "A slot with name \"" << preciceData.slotName << "\" does not exist. Available slots: " << slotNames
          << ". Use `None` or the empty string to specify the slot by its number (referenced slots will be in the order of "
          << "the items under \"preciceData\").\n Using slot No " << slotNo;
      }
    }

    // parse the mesh
    std::string preciceMeshName = currentPreciceData.getOptionString("preciceMeshName", "");
    int preciceMeshId = preciceSolverInterface_->getMeshID(preciceMeshName);

    LOG(DEBUG) << "parse item " << listIndex << " of \"preciceData\", isGeometryField: " << preciceData.isGeometryField 
      << ", opendihuMeshName: " << preciceData.opendihuMeshName << ", preciceMeshName: " << preciceMeshName 
      << ", preciceMeshId: " << preciceMeshId;

    // find the mesh among the already parsed precice meshes
    typename std::vector<std::shared_ptr<PreciceMesh>>::iterator iter
      = std::find_if(preciceMeshes_.begin(), preciceMeshes_.end(), [&preciceMeshId](std::shared_ptr<PreciceMesh> preciceMesh)
    {
      return preciceMesh->preciceMeshId == preciceMeshId;
    });

    // if the mesh is not in preciceMeshes_, create it and add it to preciceMeshes_
    if (iter == preciceMeshes_.end())
    {
      LOG(DEBUG) << "mesh was not yet initialized.";
      
      // create new precice mesh object
      std::shared_ptr<PreciceMesh> preciceMesh = std::make_shared<PreciceMesh>();
      preciceMesh->preciceMeshName = preciceMeshName;

      std::vector<Vec3> geometryValues;
      int nArrayItems = 1;

      // get the function space that is used to initialize the mapping
      std::shared_ptr<FunctionSpace> functionSpace = nullptr;
      if (preciceData.opendihuMeshName != "")
      {
        functionSpace = DihuContext::meshManager()->functionSpace<FunctionSpace>(preciceData.opendihuMeshName);
        LOG(DEBUG) << "Using opendihu mesh with name \"" << preciceData.opendihuMeshName << "\" to initialize precice mapping.";

        // get the node positions of the opendihu mesh to initialize the mapping
        preciceMesh->nNodesLocal = functionSpace->geometryField().nDofsLocalWithoutGhosts();
        functionSpace->geometryField().getValuesWithoutGhosts(geometryValues);
        
        // store opendihu mesh name
        preciceMesh->opendihuMeshName = functionSpace->meshName();
        
        // the following cannot happen
        if (preciceMesh->opendihuMeshName != preciceData.opendihuMeshName)
          LOG(FATAL) << "Initializing precice mesh from settings opendihuMeshName=\"" << preciceData.opendihuMeshName 
            << "\" and the resulting mesh is " << preciceMesh->opendihuMeshName;
      }
      else
      {
        LOG(DEBUG) << "get mesh from mesh partition, slot No " << preciceData.slotNo;
        LOG(DEBUG) << "slot connector data type: " << StringUtility::demangle(typeid(SlotConnectorDataType).name());
        
        // get the mesh partition
        std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBase
          = SlotConnectorDataHelper<SlotConnectorDataType>::getMeshPartitionBase(slotConnectorData, preciceData.slotNo, 0);

        if (!meshPartitionBase)
        {
          LOG(FATAL) << "Could not get mesh for slot No " << preciceData.slotNo;
        }
        else
          LOG(DEBUG) << "got mesh partition for slot No " << preciceData.slotNo;

        // get opendihu mesh name
        preciceMesh->opendihuMeshName = SlotConnectorDataHelper<SlotConnectorDataType>::getMeshName(slotConnectorData, preciceData.slotNo);
        
        int nDofsLocalWithoutGhosts = meshPartitionBase->nDofsLocalWithoutGhosts();
        nArrayItems = SlotConnectorDataHelper<SlotConnectorDataType>::nArrayItems(slotConnectorData, preciceData.slotNo);   // number of fibers if there are fibers

        preciceMesh->nNodesLocal = nDofsLocalWithoutGhosts * nArrayItems;

        // get the vector of values [0,1,...,nDofsLocalWithGhosts]
        const std::vector<PetscInt> &dofNosLocalWithGhosts = meshPartitionBase->dofNosLocal();
        std::vector<PetscInt> dofNosLocalWithoutGhosts(dofNosLocalWithGhosts.begin(), dofNosLocalWithGhosts.begin()+nDofsLocalWithoutGhosts);

        // loop over fibers if there are any
        for (int arrayIndex = 0; arrayIndex < nArrayItems; arrayIndex++)
        {
          static std::vector<Vec3> nodePositionsFiber;
          nodePositionsFiber.clear();
          SlotConnectorDataHelper<SlotConnectorDataType>::slotGetGeometryValues(slotConnectorData, preciceData.slotNo, arrayIndex,
                                                                                dofNosLocalWithoutGhosts, nodePositionsFiber);
          geometryValues.insert(geometryValues.end(), nodePositionsFiber.begin(), nodePositionsFiber.end());
        }

        LOG(DEBUG) << "collected " << geometryValues.size() << " node positions from the " << nArrayItems << " fibers";
      }

      // transform to contiguous memory layout for precice
      std::vector<double> geometryValuesContiguous(3*geometryValues.size());
      for (int entryNo = 0; entryNo < geometryValues.size(); entryNo++)
        for (int componentNo = 0; componentNo < 3; componentNo++)
          geometryValuesContiguous[3*entryNo + componentNo] = geometryValues[entryNo][componentNo];

      if (geometryValuesContiguous.size() != preciceMesh->nNodesLocal*3)
        LOG(FATAL) << "size mismatch: " << geometryValuesContiguous.size() << "!=" << preciceMesh->nNodesLocal*3;

      // get precice id from precice config
      preciceMesh->preciceMeshId = preciceSolverInterface_->getMeshID(preciceMesh->preciceMeshName);

      // resize buffer for vertex ids
      preciceMesh->preciceVertexIds.resize(preciceMesh->nNodesLocal);

      // give the node positions to precice and get the vertex ids
      // void precice::SolverInterface::setMeshVertices(int meshID, int size, const double *positions, int *ids)
      preciceSolverInterface_->setMeshVertices(preciceMesh->preciceMeshId, preciceMesh->nNodesLocal,
                                               geometryValuesContiguous.data(), preciceMesh->preciceVertexIds.data());

      // store newly created mesh in preciceData
      preciceData.preciceMesh = preciceMesh;

      // store the precice mesh to the vector of meshes
      preciceMeshes_.push_back(preciceMesh);

      // output message about mesh names
      std::stringstream message;
      message << "Initialized precice mesh \"" << preciceMesh->preciceMeshName 
        << "\" from opendihu mesh \"" << preciceMesh->opendihuMeshName << "\" with " << preciceMesh->nNodesLocal << " local nodes";
      if (nArrayItems > 1)
        message << " (and other corresponding meshes, in total " << nArrayItems << " fibers/compartments)";
      LOG(INFO) << message.str() << ".";
    }
    else
    {
      LOG(DEBUG) << "Use existing precice mesh with preciceMeshId " << preciceMeshId;
      preciceData.preciceMesh = *iter;
    }

    // parse mode
    std::string mode = currentPreciceData.getOptionString("mode", "");
    if (mode == "read")
    {
      preciceData.ioType = PreciceData::ioRead;
    }
    else if (mode == "write")
    {
      preciceData.ioType = PreciceData::ioWrite;
    }
    else
    {
      LOG(FATAL) << currentPreciceData << "[\"mode\"] is \"" << mode << "\", "
        << "possible values are: \"read\", \"write\".";
    }

    // parse variable name
    preciceData.preciceDataName = currentPreciceData.getOptionString("preciceDataName", "variable");

    // get precice data id
    preciceData.preciceDataId = preciceSolverInterface_->getDataID(
      preciceData.preciceDataName, preciceData.preciceMesh->preciceMeshId);

    // increment slotNo, this value is used for slots that are not identified by slotName in the config
    if (!preciceData.isGeometryField)
      slotNo++;

    // check if mesh size matches the specified slot
    // get the mesh partition
    std::shared_ptr<Partition::MeshPartitionBase> meshPartitionBase
      = SlotConnectorDataHelper<SlotConnectorDataType>::getMeshPartitionBase(slotConnectorData, preciceData.slotNo, 0);

    int nDofsLocalWithoutGhosts = meshPartitionBase->nDofsLocalWithoutGhosts();
    int nArrayItems = SlotConnectorDataHelper<SlotConnectorDataType>::nArrayItems(slotConnectorData, preciceData.slotNo);   // number of fibers if there are fibers

    if (nDofsLocalWithoutGhosts*nArrayItems != preciceData.preciceMesh->nNodesLocal)
    {
      LOG(DEBUG) << ", all available slots: " << SlotConnectorDataHelper<SlotConnectorDataType>::getString(slotConnectorData);
      LOG(FATAL) << currentPreciceData << ": Mesh does not match slot in PreciceAdapterVolumeCoupling.\n\n " 
        << "The " << (listIndex+1) << (listIndex==0? "st" : (listIndex==1? "nd" : (listIndex==2? "rd" : "th")))
        << " list item under \"preciceData\" uses slot " << preciceData.slotNo << " (slotName \"" << preciceData.slotName << "\")"
        << " and preciceMesh \"" << preciceData.preciceMesh->preciceMeshName << "\".\nThe slot has " 
        << nDofsLocalWithoutGhosts << " dofs * " << nArrayItems << " array items (fibers) = " 
        << nDofsLocalWithoutGhosts*nArrayItems << " dofs in total. The mesh was initialized from opendihu mesh \"" 
        << preciceData.preciceMesh->opendihuMeshName << "\" with " 
        << preciceData.preciceMesh->nNodesLocal << " nodes. \n" 
        << "(isGeometryField: " << std::boolalpha << preciceData.isGeometryField << ")\n"
        << "You can do the following:\n" 
        << "- Rename the preciceMesh to a unique name and update the precice config xml file accordingly.\n"
        << "- Check the solver structure file to find out the meshes that are associated with the slots.\n"
        << "  Check that the correct mesh is specified and that the precice mesh was not initialized earlier with a different opendihu mesh.\n"
        << "  Every precice mesh is initialized the first time in appears under \"preciceData\", sometimes reordering the entries can help.\n"
        << "  If the desired mesh is not available at the current solver, maybe insert a MapDofs class.";
    }
    
    LOG(INFO) << "Precice data \"" << preciceData.preciceDataName << "\" maps to " 
      << (preciceData.isGeometryField? "the geometry field of " : "") << "slot " << preciceData.slotNo 
      << " (\"" << preciceData.slotName << "\") and uses precice mesh \"" << preciceData.preciceMesh->preciceMeshName 
      << "\", which is opendihu mesh \"" << preciceData.preciceMesh->opendihuMeshName << "\" with " 
      << preciceData.preciceMesh->nNodesLocal << " local nodes.";

    // store preciceData to vector
    preciceData_.push_back(preciceData);
  }
}
#endif

}  // namespace
