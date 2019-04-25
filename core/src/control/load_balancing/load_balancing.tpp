#include "control/load_balancing/load_balancing.h"

#include "spatial_discretization/finite_element_method/finite_element_method.h"

#include <omp.h>
#include <sstream>

namespace Control
{

template<typename CellMLAdapter, typename DiffusionTimeStepping>
LoadBalancing<OperatorSplitting::Strang<TimeSteppingScheme::HeunAdaptiv<CellMLAdapter>,DiffusionTimeStepping>>::
LoadBalancing(DihuContext context) :
  LoadBalancingBase<OperatorSplitting::Strang<TimeSteppingScheme::HeunAdaptiv<CellMLAdapter>,DiffusionTimeStepping>>(context)
{
}

template<typename CellMLAdapter, typename DiffusionTimeStepping>
void LoadBalancing<OperatorSplitting::Strang<TimeSteppingScheme::HeunAdaptiv<CellMLAdapter>,DiffusionTimeStepping>>::
rebalance()
{
  LOG(TRACE) << "rebalance";

  // MPI barrier, to wait for all ranks to reach this code
  std::shared_ptr<Partition::RankSubset> rankSubsetGlobal = this->context_.partitionManager()->rankSubsetForCollectiveOperations();
  MPI_Barrier(rankSubsetGlobal->mpiCommunicator());

  // get information about finite element method object
  typedef typename DiffusionTimeStepping::DiscretizableInTime_Type FiniteElementMethodType;
  typedef typename FiniteElementMethodType::FunctionSpace FiberFunctionSpaceType;
  typedef typename DiffusionTimeStepping::Data::FieldVariableType DiffusionFieldVariableType;

  DiffusionTimeStepping &timeSteppingDiffusion = this->timeSteppingScheme_.timeStepping2();

  FiniteElementMethodType &finiteElementMethod = timeSteppingDiffusion.discretizableInTime();

  std::shared_ptr<Partition::RankSubset> rankSubsetFiber = finiteElementMethod.functionSpace()->meshPartition()->rankSubset();
  LOG(DEBUG) << "rankSubsetGlobal: " << *rankSubsetGlobal;
  LOG(DEBUG) << "rankSubsetFiber: " << *rankSubsetFiber;

  TimeSteppingScheme::HeunAdaptiv<CellMLAdapter> &timeSteppingHeun = this->timeSteppingScheme_.timeStepping1();

  // check if rebalancing is required, using timeSteppingHeun (TODO)
  // communicate between all ranks of the same fiber with rankSubsetFiber
  // communicate between all ranks with rankSubsetGlobal

  // if no rebalancing should be done, return

  // get information about local domain
  int nNodesLocalWithoutGhosts = finiteElementMethod.functionSpace()->meshPartition()->nNodesLocalWithoutGhosts();
  int nNodesGlobal = finiteElementMethod.functionSpace()->meshPartition()->nNodesGlobal();
  int nElementsLocal = finiteElementMethod.functionSpace()->meshPartition()->nElementsLocal();
  int nElementsGlobal = finiteElementMethod.functionSpace()->meshPartition()->nElementsGlobal();

  LOG(DEBUG) << "finiteElementMethod, nNodes: local: " << nNodesLocalWithoutGhosts << ", global: " << nNodesGlobal
    << ", nElements: local: " << nElementsLocal << ", global: " << nElementsGlobal;

  std::shared_ptr<DiffusionFieldVariableType> diffusionSolution = timeSteppingDiffusion.data().solution();
  LOG(DEBUG) << "finiteElementMethod field variable: " << *diffusionSolution;

  // get the components of the splitting
  CellMLAdapter &cellMLAdapter = timeSteppingHeun.discretizableInTime();

  // get information about cellml adapter
  int nInstances, nIntermediates, nParameters;
  cellMLAdapter.getNumbers(nInstances, nIntermediates, nParameters);
  LOG(DEBUG) << "cellMLAdapter has " << nInstances << " instances";

  typedef typename TimeSteppingScheme::HeunAdaptiv<CellMLAdapter>::Data::FieldVariableType CellMLFieldVariableType;
  std::shared_ptr<CellMLFieldVariableType> cellMLSolution = timeSteppingHeun.data().solution();
  LOG(DEBUG) << "cellML field variable: " << *cellMLSolution;

  // save all local values

  // get geometry field values of old mesh
  std::vector<Vec3> geometryFieldValues;
  finiteElementMethod.functionSpace()->geometryField().getValuesWithoutGhosts(geometryFieldValues);
  LOG(DEBUG) << "finiteElementMethod has geometryFieldValues: " << geometryFieldValues;

  std::vector<double> diffusionValues;
  diffusionSolution->getValuesWithoutGhosts(diffusionValues);

  const int nCellMLComponents = CellMLAdapter::nComponents();
  std::array<std::vector<double>,nCellMLComponents> cellmlValues;
  cellMLSolution->getValuesWithoutGhosts(cellmlValues);

  // determine new size of the local partition (TODO)
  int nElementsLocalNew = nElementsLocal;                      // number of elements on own rank after rebalancing

  std::array<int,1> nElementsPerDimensionLocal = {nElementsLocalNew};
  std::array<global_no_t,1> nElementsPerDimensionGlobal;
  std::array<int,1> nRanks = {rankSubsetFiber->size()};

  // create meshPartition for function space with the currently used ranks
  this->context_.partitionManager()->setRankSubsetForNextCreatedPartitioning(rankSubsetFiber);

  std::shared_ptr<Partition::MeshPartition<FiberFunctionSpaceType>> meshPartition
    = this->context_.partitionManager()->template createPartitioningStructuredLocal<FiberFunctionSpaceType>(
        nElementsPerDimensionGlobal, nElementsPerDimensionLocal, nRanks);

  int nNodesLocalWithoutGhostsNew = meshPartition->nNodesLocalWithoutGhosts();

  // communicate node positions (TODO)
  std::vector<Vec3> nodePositionsWithoutGhosts;
  nodePositionsWithoutGhosts = geometryFieldValues;
  nodePositionsWithoutGhosts.resize(nNodesLocalWithoutGhostsNew);

  // communicate new values for the diffusion problem (1 component per node) and the cellML problem (nCellMLComponents components per node) (TODO)
  std::vector<double> diffusionValuesNew = diffusionValues;
  std::array<std::vector<double>,nCellMLComponents> cellmlValuesNew = cellmlValues;

  LOG(DEBUG) << "create meshPartition, nElementsPerDimensionGlobal: " << nElementsPerDimensionGlobal;
  if (nElementsPerDimensionGlobal[0] != nElementsGlobal)
  {
    LOG(FATAL) << "number of global elements changed during rebalancing: old: " << nElementsGlobal << ", new: " << nElementsPerDimensionGlobal[0];
  }

  LOG(DEBUG) << "rebalancing number of elements: " << nElementsLocal << " -> " << nElementsLocalNew
    << ", number of nodes: " << nNodesLocalWithoutGhosts << " -> " << nNodesLocalWithoutGhostsNew;

  // create function space with updated mesh partition
  std::stringstream meshName;
  meshName << finiteElementMethod.functionSpace()->meshName() << "_";

  this->context_.partitionManager()->setRankSubsetForNextCreatedPartitioning(rankSubsetFiber);

  std::shared_ptr<FiberFunctionSpaceType> functionSpaceNew = this->context_.meshManager()->template createFunctionSpaceWithGivenMeshPartition<FiberFunctionSpaceType>(
    meshName.str(), meshPartition, nodePositionsWithoutGhosts, nElementsPerDimensionLocal, nRanks);

  finiteElementMethod = std::move(SpatialDiscretization::FiniteElementMethod<
    typename FiberFunctionSpaceType::Mesh,
    typename FiberFunctionSpaceType::BasisFunction,
    typename FiniteElementMethodType::Quadrature,
    typename FiniteElementMethodType::Term
  >(timeSteppingDiffusion.data().context(), functionSpaceNew));

  // create new CellMLAdapter with new function space
  cellMLAdapter = std::move(CellMLAdapter(cellMLAdapter, functionSpaceNew));

  // recreate data structures for cellML
  timeSteppingHeun.reset();
  timeSteppingHeun.initialize();    // retrieves function space from cellMLAdapter

  // set all local data
  timeSteppingHeun.data().solution()->setValuesWithoutGhosts(cellmlValuesNew);

  // recreate data structures for diffusion part
  timeSteppingDiffusion.reset();
  timeSteppingDiffusion.initialize();   // retrieves function space from finiteElementMethod

  // set all local data
  timeSteppingDiffusion.data().solution()->setValuesWithoutGhosts(diffusionValuesNew);
}

} // namespace
