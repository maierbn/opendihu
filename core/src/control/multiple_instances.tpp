#include "control/multiple_instances.h"

#ifdef HAVE_PAT
#include <pat_api.h>    // perftools, only available on hazel hen
#endif

#include <omp.h>
#include <sstream>

#include "data_management/multiple_instances.h"
#include "partition/partition_manager.h"
#include "utility/mpi_utility.h"
#include "control/performance_measurement.h"

namespace Control
{

template<class TimeSteppingScheme>
MultipleInstances<TimeSteppingScheme>::
MultipleInstances(DihuContext context) :
  context_(context["MultipleInstances"]), specificSettings_(context_.getPythonConfig()), data_(context_)
{
#ifdef HAVE_PAT
  std::string label = "initialization"
  PAT_region_begin(0, label.c_str());
#endif

  outputWriterManager_.initialize(context_, specificSettings_);
  
  //VLOG(1) << "MultipleInstances constructor, settings: " << specificSettings_;
  
  // extract the number of instances
  nInstances_ = specificSettings_.getOptionInt("nInstances", 1, PythonUtility::Positive);
   
  // parse all instance configs 
  std::vector<PythonConfig> instanceConfigs;
  
  // get the config for the first InstancesDataset instance from the list
  PyObject *instanceConfig = specificSettings_.getOptionListBegin<PyObject *>("instances");

  int i = 0;
  for(;
      !specificSettings_.getOptionListEnd("instances") && i < nInstances_;
      specificSettings_.template getOptionListNext<PyObject *>("instances", instanceConfig), i++)
  {
    instanceConfigs.push_back(PythonConfig(specificSettings_, "instances", instanceConfig));
    VLOG(3) << "i = " << i << ", instanceConfig = " << instanceConfig;
  }
    
  if (i < nInstances_)
  {
    LOG(ERROR) << "Could only create " << i << " instances from the given instances config python list, but nInstances = " << nInstances_;
    nInstances_ = i;
  }
    
  if (!specificSettings_.getOptionListEnd("instances"))
  {
    LOG(ERROR) << "Only " << nInstances_ << " instances were created, but more configurations are given.";
  }
  
  VLOG(1) << "MultipleInstances constructor, create Partitioning for " << nInstances_ << " instances";
  
  //MPIUtility::gdbParallelDebuggingBarrier();

  // determine all ranks of all computed instances
  std::set<int> ranksAllComputedInstances;
  nInstancesComputedGlobally_ = 0;
  std::vector<std::tuple<std::shared_ptr<Partition::RankSubset>, bool, PythonConfig>> rankSubsets(nInstances_);  // <rankSubset, computeOnThisRank, instanceConfig>

  int ownRankNoWorldCommunicator = this->context_.partitionManager()->rankNoCommWorld();
  int nRanksCommWorld = this->context_.partitionManager()->nRanksCommWorld();

  // parse the rank lists for all instances
  for (int instanceConfigNo = 0; instanceConfigNo < nInstances_; instanceConfigNo++)
  {
    PythonConfig instanceConfig = instanceConfigs[instanceConfigNo];
    std::get<2>(rankSubsets[instanceConfigNo]) = instanceConfig;
   
    // extract ranks for this instance
    if (!instanceConfig.hasKey("ranks"))
    {
      LOG(ERROR) << "Instance " << instanceConfigs << " has no \"ranks\" settings.";

      std::get<0>(rankSubsets[instanceConfigNo]) = nullptr;
      std::get<1>(rankSubsets[instanceConfigNo]) = false;
      continue;
    }
    else 
    {
      // extract rank list
      std::vector<int> ranks;
      instanceConfig.getOptionVector("ranks", ranks);
      
      VLOG(2) << "instance " << instanceConfigNo << " on ranks: " << ranks;

      // check if own rank is part of ranks list
      bool computeOnThisRank = false;

      bool computeSomewhere = false;
      for (int rank : ranks)
      {
        if (rank < nRanksCommWorld)
        {
          ranksAllComputedInstances.insert(rank);
          computeSomewhere = true;
        }
        if (rank == ownRankNoWorldCommunicator)
        {
          computeOnThisRank = true;
        }
      }

      if (computeSomewhere)
      {
        nInstancesComputedGlobally_++;
      }

      VLOG(2) << "compute on this rank: " << std::boolalpha << computeOnThisRank;

      // create rank subset
      std::shared_ptr<Partition::RankSubset> rankSubset = std::make_shared<Partition::RankSubset>(ranks.begin(), ranks.end());

      std::get<0>(rankSubsets[instanceConfigNo]) = rankSubset;
      std::get<1>(rankSubsets[instanceConfigNo]) = computeOnThisRank;
    }
  }

  // create the rank list with all computed instances
  rankSubsetAllComputedInstances_ = std::make_shared<Partition::RankSubset>(ranksAllComputedInstances.begin(), ranksAllComputedInstances.end());

  VLOG(1) << "rankSubsetAllComputedInstances: " << *rankSubsetAllComputedInstances_;

  // store the rank subset of all instances to partition manager, such that it can be retrived when the instances are generated
  this->context_.partitionManager()->setRankSubsetForCollectiveOperations(rankSubsetAllComputedInstances_);

  // log the number of instances that are computed by all ranks
  PerformanceMeasurement::setParameter("nInstancesComputedGlobally", nInstancesComputedGlobally_);

  // create all instances that are computed on the own rank
  for (int instanceConfigNo = 0; instanceConfigNo < nInstances_; instanceConfigNo++)
  {
    std::shared_ptr<Partition::RankSubset> rankSubset = std::get<0>(rankSubsets[instanceConfigNo]);
    bool computeOnThisRank = std::get<1>(rankSubsets[instanceConfigNo]);
    PythonConfig instanceConfig = std::get<2>(rankSubsets[instanceConfigNo]);

    if (!computeOnThisRank)
    {
      continue;
    }

    // store the rank subset containing only the own rank for the mesh of the current instance
    this->context_.partitionManager()->setRankSubsetForNextCreatedMesh(rankSubset);

    VLOG(1) << "create sub context for instance no " << instanceConfigNo << ", rankSubset: " << *rankSubset;
    instancesLocal_.emplace_back(context_.createSubContext(instanceConfig));
  }

  nInstancesLocal_ = instancesLocal_.size();
}

template<class TimeSteppingScheme>
void MultipleInstances<TimeSteppingScheme>::
advanceTimeSpan()
{
  // This method advances the simulation by the specified time span. It will be needed when this MultipleInstances object is part of a parent control element, like a coupling to 3D model.
  for (int i = 0; i < nInstancesLocal_; i++)
  {
    instancesLocal_[i].advanceTimeSpan();
  }
}

template<class TimeSteppingScheme>
void MultipleInstances<TimeSteppingScheme>::
setTimeSpan(double startTime, double endTime)
{
  for (int i = 0; i < nInstancesLocal_; i++)
  {
    instancesLocal_[i].setTimeSpan(startTime, endTime);
  }
}

template<class TimeSteppingScheme>
void MultipleInstances<TimeSteppingScheme>::
initialize()
{
  LOG(TRACE) << "MultipleInstances::initialize()";
  for (int i = 0; i < nInstancesLocal_; i++)
  {
    LOG(DEBUG) << "instance " << i << " initialize";
    instancesLocal_[i].initialize();
  }
  
  data_.setInstancesData(instancesLocal_);

#ifdef HAVE_PAT
  PAT_region_end(0);    // end region "initialization", id 0
#endif
}

template<class TimeSteppingScheme>
void MultipleInstances<TimeSteppingScheme>::
run()
{
  initialize();
 
  LOG(INFO) << "MultipleInstances: " << nInstancesComputedGlobally_ << " instance" << (nInstancesComputedGlobally_ != 1? "s" : "")
    << " to be computed in total.";

#ifdef HAVE_PAT
  std::string label = "computation"
  PAT_region_begin(1, label.c_str());
#endif

  //#pragma omp parallel for // does not work with the python interpreter
  for (int i = 0; i < nInstancesLocal_; i++)
  {
    if (omp_get_thread_num() == 0)
    {
      std::stringstream msg;
      msg << omp_get_thread_num() << ": running " << nInstancesLocal_ << " instances with " << omp_get_num_threads() << " OpenMP threads";
      LOG(DEBUG) << msg.str();
    }
    
    //instancesLocal_[i].reset();
    instancesLocal_[i].run();
  }
  
#ifdef HAVE_PAT
  PAT_region_end(1);    // end region "computation", id 1
#endif

  this->outputWriterManager_.writeOutput(this->data_);
}

template<class TimeSteppingScheme>
bool MultipleInstances<TimeSteppingScheme>::
knowsMeshType()
{
  // This is a dummy method that is currently not used, it is only important if we want to map between multiple data sets.
  assert(nInstances_ > 0);
  return instancesLocal_[0].knowsMeshType();
}
/*
template<class TimeSteppingScheme>
Vec &MultipleInstances<TimeSteppingScheme>::
solution()
{
  assert(nInstances_ > 0);
  return instancesLocal_[0].solution();
}*/

template<class TimeSteppingScheme>
void MultipleInstances<TimeSteppingScheme>::
reset()
{
  for (int i = 0; i < nInstancesLocal_; i++)
  {
    reset();
  }
}

template<class TimeSteppingScheme>
typename MultipleInstances<TimeSteppingScheme>::TransferableSolutionDataType MultipleInstances<TimeSteppingScheme>::
getSolutionForTransferInOperatorSplitting()
{
  std::vector<typename TimeSteppingScheme::TransferableSolutionDataType> output(nInstancesLocal_);

  for (int i = 0; i < nInstancesLocal_; i++)
  {
    output[i] = instancesLocal_[i].getSolutionForTransferInOperatorSplitting();
  }
  return output;
}

};
