#include "control/precice/partitioned_fibers.h"

#include <sstream>

namespace PreciceAdapter
{

template<class NestedSolver>
PartitionedFibers<NestedSolver>::
PartitionedFibers(DihuContext context) :
  Runnable(),
  context_(context["PartitionedFibers"]), nestedSolver_(this->context_), initialized_(false)
{
  // get python settings object from context
  this->specificSettings_ = this->context_.getPythonConfig();

  // parse options
  timeStepOutputInterval_ = specificSettings_.getOptionInt("timeStepOutputInterval", 100, PythonUtility::Positive);
}

template<class NestedSolver>
void PartitionedFibers<NestedSolver>::
initialize()
{
#ifdef HAVE_PRECICE

  // make sure that we initialize only once, in the next call, initialized_ is true
  if (initialized_)
    return;

  // initialize() will be called before the simulation starts.

  // add this solver to the solvers diagram, which is a SVG file that will be created at the end of the simulation.
  DihuContext::solverStructureVisualizer()->addSolver("PreciceAdapter::PartitionedFibers", true);   // hasInternalConnectionToFirstNestedSolver=true (the last argument) means output connector data is shared with the first subsolver

  // indicate in solverStructureVisualizer that now a child solver will be initialized
  DihuContext::solverStructureVisualizer()->beginChild();

  // call initialize of the nested solver
  nestedSolver_.initialize();

  // indicate in solverStructureVisualizer that the child solver initialization is done
  DihuContext::solverStructureVisualizer()->endChild();

  // set the outputConnectorData for the solverStructureVisualizer to appear in the solver diagram
  DihuContext::solverStructureVisualizer()->setOutputConnectorData(getOutputConnectorData());

  // get the splitting time step width which is the time step width to use here (not dt_3D which is the time window of precice)
  normalTimeStepWidth_ = nestedSolver_.nestedSolvers().instancesLocal()[0].timeStepWidth();
  currentTimeStepWidth_ = normalTimeStepWidth_;

  // check that the output connector is correct

  const int nStates = NestedSolver::CellmlAdapterType::nStates();

  std::shared_ptr<std::vector<
    std::shared_ptr<std::vector<
      std::shared_ptr<::Data::OutputConnectorData<typename NestedSolver::FunctionSpace, nStates, 1> >
    >>
  >> data = nestedSolver_.getOutputConnectorData();

  // count number of fibers
  int nFibers = 0;
  // loop over fibers
  for (int i = 0; i < data->size(); i++)
  {
    for (int j = 0; j < data->at(i)->size(); j++, nFibers++)
    {
      // get data of single fiber
      std::shared_ptr<::Data::OutputConnectorData<typename NestedSolver::FunctionSpace, nStates, 1> > fiberData
        = data->at(i)->at(j);

      // get function space of this fiber
      assert(!fiberData->variable1.empty());
      std::shared_ptr<::FunctionSpace::FunctionSpace<::Mesh::StructuredDeformableOfDimension<1>, ::BasisFunction::LagrangeOfOrder<1>>> functionSpace
        = fiberData->variable1[0].values->functionSpace();

      LOG(DEBUG) << "output connection slots:";
      LOG(DEBUG) << "  " << fiberData->variable1.size() << " in variable1:";
      for (int k = 0; k < fiberData->variable1.size(); k++)
      {
        LOG(DEBUG) << "  " << k << ": " << fiberData->variable1[k].values->name();
      }
      LOG(DEBUG) << "  " << fiberData->variable2.size() << " in variable2:";
      for (int k = 0; k < fiberData->variable2.size(); k++)
      {
        LOG(DEBUG) << "  " << k << ": " << fiberData->variable2[k].values->name();
      }
    }
  }

  // initialize precice
  const std::string solverName = "PartitionedFibers";
  const std::string configFileName = this->specificSettings_.getOptionString("preciceConfigFilename", "../precice-config.xml");

  std::shared_ptr<::FunctionSpace::FunctionSpace<::Mesh::StructuredDeformableOfDimension<1>, ::BasisFunction::LagrangeOfOrder<1>>> functionSpace = nestedSolver_.data().functionSpace();

  int rankNo = functionSpace->meshPartition()->rankSubset()->ownRankNo();
  int nRanks = functionSpace->meshPartition()->rankSubset()->size();

  preciceSolverInterface_ = std::make_unique<precice::SolverInterface>(solverName, configFileName, rankNo, nRanks);

  // store the node positions to precice
  std::string meshName = "PartitionedFibersMesh";
  preciceMeshId_ = preciceSolverInterface_->getMeshID(meshName);

  LOG(DEBUG) << "setMeshVertex to precice, one by one";

  int nDofsLocalWithoutGhosts = functionSpace->nDofsLocalWithoutGhosts();
  int nNodesTotal = nDofsLocalWithoutGhosts * nFibers;

  std::vector<double> geometryValuesPrecice(3*nNodesTotal);
  preciceVertexIds_.resize(nNodesTotal);

  // loop over fibers
  int fiberNo = 0;
  for (int i = 0; i < data->size(); i++)
  {
    for (int j = 0; j < data->at(i)->size(); j++, fiberNo++)
    {
      // get data of single fiber
      std::shared_ptr<::Data::OutputConnectorData<typename NestedSolver::FunctionSpace, nStates, 1> > fiberData
        = data->at(i)->at(j);

      // get function space of this fiber
      assert(!fiberData->variable1.empty());
      std::shared_ptr<::FunctionSpace::FunctionSpace<::Mesh::StructuredDeformableOfDimension<1>, ::BasisFunction::LagrangeOfOrder<1>>> functionSpace
        = fiberData->variable1[0].values->functionSpace();

      LOG(DEBUG) << "fiberData " << i << " " << j << ": " << fiberData->variable1.size() << "," << fiberData->variable2.size() << ": " << *fiberData;

      // loop over dofs of fiber
      for (int dofNoLocal = 0; dofNoLocal < nDofsLocalWithoutGhosts; dofNoLocal++)
      {
        // set node position to precice
        Vec3 nodePosition = functionSpace->getGeometry(dofNoLocal);

        for (int componentNo = 0; componentNo < 3; componentNo++)
        {
          geometryValuesPrecice[3*(fiberNo*nDofsLocalWithoutGhosts + dofNoLocal) + componentNo] = nodePosition[componentNo];
        }
      }
    }
  }

  preciceSolverInterface_->setMeshVertices(preciceMeshId_, nNodesTotal, geometryValuesPrecice.data(), preciceVertexIds_.data());

  LOG(DEBUG) << "precice defined vertexIds: " << preciceVertexIds_;

  // initialize data ids
  preciceDataIdGeometry_  = preciceSolverInterface_->getDataID("Geometry",  preciceMeshId_);
  preciceDataIdGamma_     = preciceSolverInterface_->getDataID("Gamma",     preciceMeshId_);
  LOG(DEBUG) << "data id geometry: " << preciceDataIdGeometry_;
  LOG(DEBUG) << "data id gamma: " << preciceDataIdGamma_;
  //preciceDataIdLambda_    = preciceSolverInterface_->getDataID("Lambda",    preciceMeshId_);
  //preciceDataIdLambdaDot_ = preciceSolverInterface_->getDataID("LambdaDot", preciceMeshId_);

  //preciceSolverInterface_->setMeshVertices(preciceMeshId_, int size, double* positions, int* ids);


  maximumPreciceTimestepSize_ = preciceSolverInterface_->initialize();

  LOG(DEBUG) << "precice initialization done, dt: " << maximumPreciceTimestepSize_;

  initialized_ = true;

#else
  LOG(FATAL) << "Not compiled with preCICE!";
#endif
}

template<class NestedSolver>
void PartitionedFibers<NestedSolver>::
run()
{
#ifdef HAVE_PRECICE

  // initialize everything
  initialize();

  // data to send:
  // - gamma
  // data to receive:
  // - geometry
  // - lambda
  // - lambdaDot

  // nestedSolver_:
  // FastMonodomainSolver<Control::MultipleInstances<OperatorSplitting::Strang<Control::MultipleInstances<TimeSteppingScheme::Heun<CellmlAdapter<57, 71, FunctionSpace::FunctionSpace<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1> > > > >, Control::MultipleInstances<TimeSteppingScheme::ImplicitEuler<SpatialDiscretization::FiniteElementMethod<Mesh::StructuredDeformableOfDimension<1>, BasisFunction::LagrangeOfOrder<1>, Quadrature::Gauss<2u>, Equation::Dynamic::IsotropicDiffusion, Equation::Dynamic::IsotropicDiffusion, BasisFunction::LagrangeOfOrder<1> > > > > > >

  if (preciceSolverInterface_->isActionRequired(precice::constants::actionWriteInitialData()))
  {
    // writeData
    preciceWriteData();

    preciceSolverInterface_->markActionFulfilled(precice::constants::actionWriteInitialData());

    // initialize data in precice
    preciceSolverInterface_->initializeData();
  }

  double currentSimulationTime = 0;

  // main simulation loop of adapter
  for (int timeStepNo = 0; preciceSolverInterface_->isCouplingOngoing(); timeStepNo++)
  {
    // read data from precice
    // data to receive:
    // - geometry
    // - lambda
    // - lambdaDot
    preciceReadData();

    // compute the time step width such that it fits in the remaining time in the current time window
    currentTimeStepWidth_ = std::min(maximumPreciceTimestepSize_, normalTimeStepWidth_);

    // set the timestep width in the solver
    nestedSolver_.setTimeSpan(currentSimulationTime, currentSimulationTime+currentTimeStepWidth_);

    if (timeStepNo % this->timeStepOutputInterval_ == 0 && (this->timeStepOutputInterval_ <= 10 || timeStepNo > 0))  // show first timestep only if timeStepOutputInterval is <= 10
    {
      LOG(INFO) << "PartitionedFibers, timestep " << timeStepNo << ", t=" << currentSimulationTime
        << ", dt: " << currentTimeStepWidth_;
    }

    // call the nested solver
    nestedSolver_.advanceTimeSpan();

    // write data to precice
    // data to send:
    // - gamma
    preciceWriteData();

    currentSimulationTime += currentTimeStepWidth_;

    LOG(DEBUG) << "precice::advance(" << currentTimeStepWidth_ << ")";

    // advance precice
    maximumPreciceTimestepSize_ = preciceSolverInterface_->advance(currentTimeStepWidth_);
  }

  preciceSolverInterface_->finalize();
#endif
}

#ifdef HAVE_PRECICE

template<class NestedSolver>
void PartitionedFibers<NestedSolver>::
preciceReadData()
{
  if (!preciceSolverInterface_->isReadDataAvailable())
    return;

  LOG(DEBUG) << "read data from precice";

  // read data from precice, only if not first timestep

  // data to receive:
  // - geometry
  // - lambda
  // - lambdaDot
  const int nStates = NestedSolver::CellmlAdapterType::nStates();

  std::shared_ptr<std::vector<
    std::shared_ptr<std::vector<
      std::shared_ptr<::Data::OutputConnectorData<typename NestedSolver::FunctionSpace, nStates, 1> >
    >>
  >> data = nestedSolver_.getOutputConnectorData();

  // loop over fibers
  int fiberNo = 0;
  for (int i = 0; i < data->size(); i++)
  {
    for (int j = 0; j < data->at(i)->size(); j++, fiberNo++)
    {
      // get data for a single fiber
      std::shared_ptr<::Data::OutputConnectorData<typename NestedSolver::FunctionSpace, nStates, 1> > fiberData
        = data->at(i)->at(j);

      // get function space
      assert(!fiberData->variable1.empty());
      std::shared_ptr<::FunctionSpace::FunctionSpace<::Mesh::StructuredDeformableOfDimension<1>, ::BasisFunction::LagrangeOfOrder<1>>> functionSpace
        = fiberData->variable1[0].values->functionSpace();

      // loop over dofs of a single fiber mesh
      for (int dofNoLocal = 0; dofNoLocal < functionSpace->nDofsLocalWithoutGhosts(); dofNoLocal++)
      {
        // read geometry value from precice
        Vec3 nodePosition;
        assert(fiberNo * functionSpace->nDofsLocalWithoutGhosts() + dofNoLocal < preciceVertexIds_.size());

        int vertexIndex = preciceVertexIds_[fiberNo * functionSpace->nDofsLocalWithoutGhosts() + dofNoLocal];
        preciceSolverInterface_->readVectorData(preciceDataIdGeometry_, vertexIndex, nodePosition.data());

        // store node position
        functionSpace->geometryField().setValue(dofNoLocal, nodePosition);
#if 0
        // read scalar lambda value from precice
        double value;
        preciceSolverInterface_->readScalarData(preciceDataIdLambda_, vertexIndex, value);

        // set lambda value
        fiberData->variable1[0].setValue(dofNoLocal, value);

        // read scalar lambdaDot value from precice
        preciceSolverInterface_->readScalarData(preciceDataIdLambda_, vertexIndex, value);

        // set lambdaDot value
        assert(fiberData->variable1.size() > 1);
        fiberData->variable1[1].setValue(dofNoLocal, value);
#endif
      }   // dofNoLocal
    }   // j
  }   // i

  LOG(DEBUG) << "precice read geometry complete";
}

template<class NestedSolver>
void PartitionedFibers<NestedSolver>::
preciceWriteData()
{
  if (!preciceSolverInterface_->isWriteDataRequired(currentTimeStepWidth_))
    return;

  LOG(DEBUG) << "write gamma data to precice, one by one";

  // write data to precice
  // data to send:
  // - gamma

  const int nStates = NestedSolver::CellmlAdapterType::nStates();

  std::shared_ptr<std::vector<
    std::shared_ptr<std::vector<
      std::shared_ptr<::Data::OutputConnectorData<typename NestedSolver::FunctionSpace, nStates, 1> >
    >>
  >> data = nestedSolver_.getOutputConnectorData();

  // loop over fibers
  int fiberNo = 0;
  for (int i = 0; i < data->size(); i++)
  {
    for (int j = 0; j < data->at(i)->size(); j++, fiberNo++)
    {
      // get data for a single fiber
      std::shared_ptr<::Data::OutputConnectorData<typename NestedSolver::FunctionSpace, nStates, 1> > fiberData
        = data->at(i)->at(j);

        // get function space
      assert(!fiberData->variable1.empty());
      std::shared_ptr<::FunctionSpace::FunctionSpace<::Mesh::StructuredDeformableOfDimension<1>, ::BasisFunction::LagrangeOfOrder<1>>> functionSpace
        = fiberData->variable1[0].values->functionSpace();

      // loop over dofs of a single fiber mesh
      for (int dofNoLocal = 0; dofNoLocal < functionSpace->nDofsLocalWithoutGhosts(); dofNoLocal++)
      {
        // get gamma value
        double gamma = fiberData->variable2[0].getValue(dofNoLocal);

        // write gamma value to precice
        assert(fiberNo * functionSpace->nDofsLocalWithoutGhosts() + dofNoLocal < preciceVertexIds_.size());

        int vertexIndex = preciceVertexIds_[fiberNo * functionSpace->nDofsLocalWithoutGhosts() + dofNoLocal];
        preciceSolverInterface_->writeScalarData(preciceDataIdGamma_, vertexIndex, gamma);
      }   // dofNoLocal
    }   // j
  }   // i

  LOG(DEBUG) << "write gamma data to precice complete";
}

#endif

template<class NestedSolver>
void PartitionedFibers<NestedSolver>::
reset()
{
  nestedSolver_.reset();

  initialized_ = false;
  // "uninitialize" everything
}

template<class NestedSolver>
typename PartitionedFibers<NestedSolver>::Data &PartitionedFibers<NestedSolver>::
data()
{
  // get a reference to the data object
  return nestedSolver_.data();
}

//! get the data that will be transferred in the operator splitting to the other term of the splitting
//! the transfer is done by the output_connector_data_transfer class
template<class NestedSolver>
std::shared_ptr<typename PartitionedFibers<NestedSolver>::OutputConnectorDataType> PartitionedFibers<NestedSolver>::
getOutputConnectorData()
{
  //! This is relevant only, if this solver is part of a splitting or coupling scheme. Then this method returns the values/variables that will be
  // transferred to the other solvers. We can just reuse the values of the nestedSolver_.
  return nestedSolver_.getOutputConnectorData();
}

}  // namespace
