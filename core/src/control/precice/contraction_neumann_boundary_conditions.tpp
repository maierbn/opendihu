#include "control/precice/contraction_neumann_boundary_conditions.h"

#include "spatial_discretization/neumann_boundary_conditions/01_neumann_boundary_conditions.h"

#include <sstream>

namespace PreciceAdapter
{

template<typename NestedSolver>
ContractionNeumannBoundaryConditions<NestedSolver>::
ContractionNeumannBoundaryConditions(DihuContext context) :
  Runnable(),
  context_(context["PreciceContractionNeumannBoundaryConditions"]), nestedSolver_(this->context_), initialized_(false)
{
  // get python settings object from context
  this->specificSettings_ = this->context_.getPythonConfig();
}

template<typename NestedSolver>
void ContractionNeumannBoundaryConditions<NestedSolver>::
initialize()
{
#ifdef HAVE_PRECICE

  // make sure that we initialize only once, in the next call, initialized_ is true
  if (initialized_)
    return;

  // initialize() will be called before the simulation starts.

  // add this solver to the solvers diagram, which is an ASCII art representation that will be created at the end of the simulation.
  DihuContext::solverStructureVisualizer()->addSolver("PreciceAdapter::ContractionNeumannBoundaryConditions", true);   // hasInternalConnectionToFirstNestedSolver=true (the last argument) means output connector data is shared with the first subsolver

  // indicate in solverStructureVisualizer that now a child solver will be initialized
  DihuContext::solverStructureVisualizer()->beginChild();

  // call initialize of the nested solver
  nestedSolver_.initialize();

  // indicate in solverStructureVisualizer that the child solver initialization is done
  DihuContext::solverStructureVisualizer()->endChild();

  // initialize precice
  const std::string solverName = "TendonSolver";
  const std::string configFileName = this->specificSettings_.getOptionString("preciceConfigFilename", "../precice-config.xml");

  functionSpace_ = nestedSolver_.data().functionSpace();
  int rankNo = functionSpace_->meshPartition()->rankSubset()->ownRankNo();
  int nRanks = functionSpace_->meshPartition()->rankSubset()->size();

  preciceSolverInterface_ = std::make_unique<precice::SolverInterface>(solverName, configFileName, rankNo, nRanks);

  // store the node positions to precice
  std::string meshName = "TendonTopMesh";
  preciceMeshId_ = preciceSolverInterface_->getMeshID(meshName);


  std::vector<Vec3> geometryValues;
  functionSpace_->geometryField().getValuesWithoutGhosts(geometryValues);

  // get nodes at coupling surface
  const int nNodesX = functionSpace_->nNodesLocalWithoutGhosts(0);
  const int nNodesY = functionSpace_->nNodesLocalWithoutGhosts(1);
  nNodesSurfaceLocal_ = nNodesX * nNodesY;

  std::vector<double> geometryValuesSurfacePrecice(3*nNodesSurfaceLocal_);
  preciceVertexIds_.resize(nNodesSurfaceLocal_);

  // get the parameter if the coupling surface is at z- or at z+
  isCouplingSurfaceBottom_ = this->specificSettings_.getOptionBool("isCouplingSurfaceBottom", false);

  int nodeIndexZ = 0;
  if (isCouplingSurfaceBottom_)
    nodeIndexZ = functionSpace_->nNodesLocalWithoutGhosts(2) - 1;

  // loop over nodes
  for (int nodeIndexY = 0; nodeIndexY < nNodesY; nodeIndexY++)
  {
    for (int nodeIndexX = 0; nodeIndexX < nNodesX; nodeIndexX++)
    {
      node_no_t nodeNoLocal =
        nodeIndexZ * nNodesX * nNodesY
        + nodeIndexY * nNodesX
        + nodeIndexX;
      dof_no_t dofNoLocal = nodeNoLocal;

      int valueIndex = nodeIndexY*nNodesX + nodeIndexX;

      for (int i = 0; i < 3; i++)
      {
        geometryValuesSurfacePrecice[3*valueIndex + i] = geometryValues[dofNoLocal][i];
      }
    }
  }

  LOG(DEBUG) << "setMeshVertices to precice, " << 3*nNodesSurfaceLocal_ << " values";

  //preciceSolverInterface_->setMeshVertices(preciceMeshId_, int size, double* positions, int* ids);
  preciceSolverInterface_->setMeshVertices(preciceMeshId_, nNodesSurfaceLocal_, geometryValuesSurfacePrecice.data(), preciceVertexIds_.data());

  LOG(DEBUG) << "precice defined vertexIds: " << preciceVertexIds_;

  // initialize data ids
  preciceDataIdDisplacements_ = preciceSolverInterface_->getDataID("Displacements", preciceMeshId_);
  preciceDataIdVelocity_      = preciceSolverInterface_->getDataID("Velocity",      preciceMeshId_);
  preciceDataIdTraction_      = preciceSolverInterface_->getDataID("Traction",      preciceMeshId_);

  LOG(DEBUG) << "data id displacements: " << preciceDataIdDisplacements_;
  LOG(DEBUG) << "data id traction: " << preciceDataIdTraction_;

  maximumPreciceTimestepSize_ = preciceSolverInterface_->initialize();

  timeStepWidth_ = this->specificSettings_.getOptionDouble("timestepWidth", 0.01, PythonUtility::Positive);
  LOG(DEBUG) << "precice initialization done, dt: " << maximumPreciceTimestepSize_ << "," << timeStepWidth_;

  initialized_ = true;

#else
  LOG(FATAL) << "Not compiled with preCICE!";
#endif
}

template<typename NestedSolver>
void ContractionNeumannBoundaryConditions<NestedSolver>::
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

  //std::shared_ptr<::Data::OutputConnectorData<::FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<3>, ::BasisFunction::LagrangeOfOrder<2> >, 1, 1> > connectorData
  //  = getOutputConnectorData();

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
    // - displacements
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
void ContractionNeumannBoundaryConditions<NestedSolver>::
preciceReadData()
{
  if (!preciceSolverInterface_->isReadDataAvailable())
    return;

  LOG(DEBUG) << "read data from precice";

  // read traction values from precice
  std::vector<double> tractionValues(nNodesSurfaceLocal_*3);
  preciceSolverInterface_->readBlockVectorData(preciceDataIdTraction_, nNodesSurfaceLocal_,
                                               preciceVertexIds_.data(), tractionValues.data());
  // set traction values as neumann boundary conditions
  using ElementWithFacesType = typename SpatialDiscretization::NeumannBoundaryConditions<FunctionSpace,Quadrature::Gauss<3>,3>::ElementWithFaces;
  /*
  struct ElementWithFaces
  {
    element_no_t elementNoLocal;                                     //< the local no of the element

    Mesh::face_t face;                                               //< face on which the Neumann BC is applied
    std::vector<std::pair<dof_no_t, VecD<nComponents>>> dofVectors;  //< <surface-local dof no, value>, nComponents == FunctionSpaceType::dim() for traction boundary condition or nComponents = 1 for flux BC
    std::vector<dof_no_t> surfaceDofs;                               //< dof nos of the volume element that correspond to the face / surface. These are different from the dofs in dofsVector which are numbered for the surface only, surfaceDofs are in the numbering of the volume element.
    // note, for flux BC, dofVectors[i].second is a VecD<1>
  };*/


  std::vector<ElementWithFacesType> neumannBoundaryConditionElements;

  const int nNodesX = functionSpace_->nNodesLocalWithoutGhosts(0);
  const int nNodesY = functionSpace_->nNodesLocalWithoutGhosts(1);
  const int nNodesZ = functionSpace_->nNodesLocalWithoutGhosts(2);

  int nodeIndexZ = nNodesZ-1;
  if (isCouplingSurfaceBottom_)
    nodeIndexZ = 0;

  int nElementsX = functionSpace_->meshPartition()->nElementsLocal(0);
  int nElementsY = functionSpace_->meshPartition()->nElementsLocal(1);
  int nElementsZ = functionSpace_->meshPartition()->nElementsLocal(2);
  int elementIndexZ = nElementsZ-1;
  if (isCouplingSurfaceBottom_)
    elementIndexZ = 0;

  // loop over elements
  for (int elementIndexY = 0; elementIndexY < nElementsY; elementIndexY++)
  {
    for (int elementIndexX = 0; elementIndexX < nElementsX; elementIndexX++)
    {
      ElementWithFacesType elementWithFaces;
      element_no_t elementNoLocal = elementIndexZ*nElementsX*nElementsY + elementIndexY * nElementsX + elementIndexX;
      elementWithFaces.elementNoLocal = elementNoLocal;

      // set surface dofs
      Mesh::face_t face = Mesh::face_t::face2Plus;
      if (isCouplingSurfaceBottom_)
        face = Mesh::face_t::face2Minus;

      // get dofs indices within the numbering of the volume element that correspond to the selected face
      const int nDofsPerNode = FunctionSpace::nDofsPerNode();
      const int nSurfaceDofs = ::FunctionSpace::FunctionSpaceBaseDim<2,typename FunctionSpace::BasisFunction>::nNodesPerElement() * nDofsPerNode;
      std::array<dof_no_t,nSurfaceDofs> surfaceDofs;
      FunctionSpace::getFaceDofs(face, surfaceDofs);
      elementWithFaces.surfaceDofs.assign(surfaceDofs.begin(), surfaceDofs.end());

      // loop over the nodes of the element
      for (int indexY = 0; indexY < 3; indexY++)
      {
        for (int indexX = 0; indexX < 3; indexX++)
        {
          int elementalDofIndex = indexY * 3 + indexX;

          dof_no_t dofNoLocal = functionSpace_->getDofNo(elementNoLocal, elementalDofIndex);
          int valueIndex = dofNoLocal - nodeIndexZ*nNodesX*nNodesY;

          Vec3 traction;
          for (int i = 0; i < 3; i++)
          {
            traction[i] = tractionValues[3*valueIndex + i];
          }

          if (isCouplingSurfaceBottom_)
          {
            dof_no_t surfaceDof = elementalDofIndex;
            elementWithFaces.dofVectors.push_back(std::pair<dof_no_t,Vec3>(surfaceDof, traction));
          }
          else
          {
            dof_no_t surfaceDof = 18+elementalDofIndex;
            elementWithFaces.dofVectors.push_back(std::pair<dof_no_t,Vec3>(surfaceDof, traction));
          }
        }
      }
      neumannBoundaryConditionElements.push_back(elementWithFaces);
    }
  }

  // create new Neumann BC object
  using NeumannBoundaryConditionsType = SpatialDiscretization::NeumannBoundaryConditions<FunctionSpace,Quadrature::Gauss<3>,3>;
  std::shared_ptr<NeumannBoundaryConditionsType> neumannBoundaryConditions = std::make_shared<NeumannBoundaryConditionsType>(this->context_);
  neumannBoundaryConditions->initialize(functionSpace_, neumannBoundaryConditionElements);

  nestedSolver_.hyperelasticitySolver().updateNeumannBoundaryConditions(neumannBoundaryConditions);

  //getOutputConnectorData()->variable1[outputConnectorSlotIdGamma_].setValuesWithoutGhosts(gammaValues);

  LOG(DEBUG) << "read data from precice complete, traction values: " << tractionValues;
}

template<typename NestedSolver>
void ContractionNeumannBoundaryConditions<NestedSolver>::
preciceWriteData()
{
  if (!preciceSolverInterface_->isWriteDataRequired(timeStepWidth_))
    return;

  // write displacements data to precice

  // convert geometry values to precice data layout
  std::vector<Vec3> displacementsValues;
  nestedSolver_.data().displacements()->getValuesWithoutGhosts(displacementsValues);
  std::vector<double> displacementValuesPrecice(3*nNodesSurfaceLocal_);

  std::vector<Vec3> velocityValues;
  nestedSolver_.data().velocities()->getValuesWithoutGhosts(velocityValues);
  std::vector<double> velocityValuesPrecice(3*nNodesSurfaceLocal_);

  // loop over top nodes
  const int nNodesX = functionSpace_->nNodesLocalWithoutGhosts(0);
  const int nNodesY = functionSpace_->nNodesLocalWithoutGhosts(1);
  const int nNodesZ = functionSpace_->nNodesLocalWithoutGhosts(2);

  int nodeIndexZ = nNodesZ-1;
  if (isCouplingSurfaceBottom_)
    nodeIndexZ = 0;

  for (int nodeIndexY = 0; nodeIndexY < nNodesY; nodeIndexY++)
  {
    for (int nodeIndexX = 0; nodeIndexX < nNodesX; nodeIndexX++)
    {
      node_no_t nodeNoLocal =
        nodeIndexZ * nNodesX * nNodesY
        + nodeIndexY * nNodesX
        + nodeIndexX;

      dof_no_t dofNoLocal = nodeNoLocal;
      int valueIndex = nodeIndexY*nNodesX + nodeIndexX;

      for (int i = 0; i < 3; i++)
      {
        displacementValuesPrecice[3*valueIndex + i] = displacementsValues[dofNoLocal][i];
        velocityValuesPrecice[3*valueIndex + i] = velocityValues[dofNoLocal][i];
      }
    }
  }

  LOG(DEBUG) << "write displacements data to precice: " << displacementValuesPrecice;

  // write displacement values in precice
  preciceSolverInterface_->writeBlockVectorData(preciceDataIdDisplacements_, nNodesSurfaceLocal_,
                                                preciceVertexIds_.data(), displacementValuesPrecice.data());

  // write velocity values in precice
  preciceSolverInterface_->writeBlockVectorData(preciceDataIdVelocity_, nNodesSurfaceLocal_,
                                                preciceVertexIds_.data(), velocityValuesPrecice.data());

  LOG(DEBUG) << "write displacements data to precice complete";
}

#endif

template<typename NestedSolver>
void ContractionNeumannBoundaryConditions<NestedSolver>::
reset()
{
  nestedSolver_.reset();

  initialized_ = false;
  // "uninitialize" everything
}

template<typename NestedSolver>
typename ContractionNeumannBoundaryConditions<NestedSolver>::Data &ContractionNeumannBoundaryConditions<NestedSolver>::
data()
{
  // get a reference to the data object
  return nestedSolver_.data();
}

}  // namespace
