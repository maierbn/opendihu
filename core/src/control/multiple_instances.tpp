#include "control/multiple_instances.h"

#ifdef HAVE_PAT
#include <pat_api.h>    // perftools, only available on hazel hen
#endif

#include <omp.h>
#include <sstream>
#include <string>

#include "data_management/control/multiple_instances.h"
#include "partition/partition_manager.h"
#include "utility/mpi_utility.h"
#include "control/diagnostic_tool/performance_measurement.h"
#include "control/diagnostic_tool/solver_structure_visualizer.h"

namespace Control
{

template<typename TimeSteppingScheme>
MultipleInstances<TimeSteppingScheme>::
MultipleInstances(DihuContext context) :
  context_(context["MultipleInstances"]), specificSettings_(context_.getPythonConfig()),
  data_(context_), outputInitializeThisInstance_(false)
{
  std::vector<std::string> configKeys;
  specificSettings_.getKeys(configKeys);
  LOG(DEBUG) << "initialize outputWriterManager_, keys: " << configKeys;

  // parse log key
  if (specificSettings_.hasKey("durationLogKey") || specificSettings_.hasKey("logKey"))
  {
    if (specificSettings_.hasKey("logKey"))
    {
      this->logKey_ = specificSettings_.getOptionString("logKey", "");
    }
    else
    {
      this->logKey_ = specificSettings_.getOptionString("durationLogKey", "");
    }
  }

  outputWriterManager_.initialize(context_, specificSettings_);
  
  //LOG(DEBUG) << "MultipleInstances constructor, settings: ";
  //PythonUtility::printDict(specificSettings_.pyObject());
  
  // extract the number of instances
  nInstances_ = specificSettings_.getOptionInt("nInstances", 1, PythonUtility::Positive);

  // parse all instance configs 
  std::vector<std::shared_ptr<PythonConfig>> instanceConfigs;
  
  // get the config for the first InstancesDataset instance from the list
  PyObject *instanceConfig = specificSettings_.getOptionListBegin<PyObject *>("instances");

  int i = 0;
  for (;
      !specificSettings_.getOptionListEnd("instances") && i < nInstances_;
      specificSettings_.template getOptionListNext<PyObject *>("instances", instanceConfig), i++)
  {
    if (instanceConfig == Py_None)
    {
      instanceConfigs.push_back(nullptr);
    }
    else
    {
      instanceConfigs.push_back(std::make_shared<PythonConfig>(specificSettings_, "instances", instanceConfig, i));
    }

    VLOG(3) << "i = " << i << ", instanceConfig = " << instanceConfig;
  }
    
  if (i < nInstances_)
  {
    LOG(ERROR) << "Could only create " << i << " instances from the given instances config python list, but nInstances = " << nInstances_;
    nInstances_ = i;
  }
    
  if (!specificSettings_.getOptionListEnd("instances"))
  {
    PyObject *instancesList = specificSettings_.getOptionPyObject("instances");
    std::vector<PyObject*> vector = PythonUtility::convertFromPython<std::vector<PyObject *>>::get(instancesList);
    LOG(ERROR) << "Only " << nInstances_ << " instances were created, but more (" << vector.size() << ") configurations are given.";
  }
  
  VLOG(1) << "MultipleInstances constructor, create Partitioning for " << nInstances_ << " instances";
  
  //MPIUtility::gdbParallelDebuggingBarrier();

  // determine all ranks of all computed instances
  std::set<int> ranksAllComputedInstances;

  // if all ranks are given, parse them
  if (specificSettings_.hasKey("ranksAllComputedInstances"))
  {
    std::vector<int> ranksAllComputedInstancesVector;
    specificSettings_.getOptionVector("ranksAllComputedInstances", ranksAllComputedInstancesVector);
    ranksAllComputedInstances.insert(ranksAllComputedInstancesVector.begin(), ranksAllComputedInstancesVector.end());
  }

  nInstancesComputedGlobally_ = 0;
  std::vector<std::tuple<std::shared_ptr<Partition::RankSubset>, bool, std::shared_ptr<PythonConfig>>> instanceData(nInstances_);  // <rankSubset, computeOnThisRank, instanceConfig>

  int ownRankNo = this->context_.ownRankNo();  // this may not be from MPI_COMM_WORLD but the context's communicator
  int nRanksThisContext = this->context_.rankSubset()->size();
  std::vector<std::shared_ptr<Partition::RankSubset>> rankSubsets;

  // parse the rank lists for all instances
  for (int instanceConfigNo = 0; instanceConfigNo < nInstances_; instanceConfigNo++)
  {
    std::shared_ptr<PythonConfig> instanceConfig = instanceConfigs[instanceConfigNo];
    std::get<2>(instanceData[instanceConfigNo]) = instanceConfig;
   
    std::shared_ptr<Partition::RankSubset> rankSubset = nullptr;
    std::vector<int> ranks;
    bool computeOnThisRank = false;

    // extract ranks for this instance
    if (!instanceConfig)
    {
      // do nothing
    }
    else if (!instanceConfig->hasKey("ranks"))
    {
      LOG(ERROR) << "Instance " << instanceConfigs << " has no \"ranks\" settings.";
    }
    else 
    {
      // extract rank list from config
      instanceConfig->getOptionVector("ranks", ranks);
      std::set<int> rankSet(ranks.begin(), ranks.end());
      LOG(DEBUG) << "instance no. " << instanceConfigNo << " has ranks: " << ranks << ", " << rankSet;

      VLOG(2) << "instance " << instanceConfigNo << " on ranks: " << ranks;

      // check if own rank is part of ranks list
      bool computeSomewhere = false;
      for (int rank : rankSet)
      {
        if (rank < nRanksThisContext)
        {
          ranksAllComputedInstances.insert(rank);
          computeSomewhere = true;
        }
        if (rank == ownRankNo)
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
      if (computeOnThisRank)
      {
        LOG(DEBUG) << "instance " << instanceConfigNo << ", compute on this rank";
      }

      // check if a matching rank subset already exists that can be reused
      for (int i = 0; i < rankSubsets.size(); i++)
      {
        if (rankSubsets[i]->equals(rankSet))
        {
          rankSubset = rankSubsets[i];
          LOG(DEBUG) << "reuse rank subset of instance " << i << ", ranks: " << ranks;
          break;
        }
      }
    }

    if (!rankSubset)
    {
      LOG(DEBUG) << "create new rank subset from ranks " << ranks << " in context " << *this->context_.rankSubset();
      // The rank subsets have to be created collectively by all ranks in the current context, even if they will not be part of the new communicator!
      rankSubset = std::make_shared<Partition::RankSubset>(ranks.begin(), ranks.end(), this->context_.rankSubset());
    }
    rankSubsets.push_back(rankSubset);

    std::get<0>(instanceData[instanceConfigNo]) = rankSubset;
    std::get<1>(instanceData[instanceConfigNo]) = computeOnThisRank;
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
    std::shared_ptr<Partition::RankSubset> rankSubset = std::get<0>(instanceData[instanceConfigNo]);
    bool computeOnThisRank = std::get<1>(instanceData[instanceConfigNo]);
    std::shared_ptr<PythonConfig> instanceConfig = std::get<2>(instanceData[instanceConfigNo]);

    if (!computeOnThisRank)
    {
      continue;
    }

    // store the rank subset containing only the own rank for the mesh of the current instance
    this->context_.partitionManager()->setRankSubsetForNextCreatedPartitioning(rankSubset);

    VLOG(1) << "create sub context for instance no " << instanceConfigNo << ", rankSubset: " << *rankSubset;
    VLOG(1) << "this rankSubset is also set via setRankSubsetForNextCreatedPartitioning";

    // add this instance to the instances that are computed locally
    instancesLocal_.emplace_back(context_.createSubContext(*instanceConfig, rankSubset));
    rankSubsetsLocal_.push_back(rankSubset);
  }

  nInstancesLocal_ = instancesLocal_.size();

  // store the number of local instances to be included in the log file
  if (this->logKey_ != "")
  {
    std::stringstream logKey;
    logKey << this->logKey_ << "_n";
    Control::PerformanceMeasurement::setParameter(logKey.str(), nInstancesLocal_);
  }

  // clear rank subset for next created partitioning
  VLOG(1) << "clear rank subset for next created partitioning";
  this->context_.partitionManager()->setRankSubsetForNextCreatedPartitioning(nullptr);
}

template<typename TimeSteppingScheme>
void MultipleInstances<TimeSteppingScheme>::
advanceTimeSpan(bool withOutputWritersEnabled)
{
  // start duration measurement
  if (this->logKey_ != "")
    Control::PerformanceMeasurement::start(this->logKey_);

  // This method advances the simulation by the specified time span. It will be needed when this MultipleInstances object is part of a parent control element, like a coupling to 3D model.
  for (int i = 0; i < nInstancesLocal_; i++)
  {
    instancesLocal_[i].advanceTimeSpan(withOutputWritersEnabled);
  }

  // stop duration measurement
  if (this->logKey_ != "")
    Control::PerformanceMeasurement::stop(this->logKey_);

  LOG(DEBUG) << "multipleInstances::advanceTimeSpan() complete, now call writeOutput, hasOutputWriters: " << this->outputWriterManager_.hasOutputWriters();

  if (nInstancesLocal_ > 0)
  {
    writeOwnOutput(instancesLocal_[0].numberTimeSteps(), instancesLocal_[0].endTime());
  }
}

template<typename TimeSteppingScheme>
void MultipleInstances<TimeSteppingScheme>::
setTimeSpan(double startTime, double endTime)
{
  for (int i = 0; i < nInstancesLocal_; i++)
  {
    instancesLocal_[i].setTimeSpan(startTime, endTime);
  }
}

template<typename TimeSteppingScheme>
void MultipleInstances<TimeSteppingScheme>::
initialize()
{
  if (this->logKey_ != "")
  {
    std::stringstream logKey;
    logKey << this->logKey_ << "_init";
    Control::PerformanceMeasurement::start(logKey.str());
  }

  LOG(TRACE) << "MultipleInstances::initialize()";

  // initialize output of progress in %, it is only output for one instance and then only for rank 0
  if (outputInitialize_)
  {
    outputInitializeThisInstance_ = true;
    outputInitialize_ = false;
    el::Loggers::removeFlag(el::LoggingFlag::NewLineForContainer);
    LOG(INFO) << "Initialize " << nInstancesComputedGlobally_ << " global instances (" << nInstancesLocal_ << " local).";
  }

  // add this solver to the solvers diagram
  DihuContext::solverStructureVisualizer()->addSolver("MultipleInstances", true);   // hasInternalConnectionToFirstNestedSolver=true (the last argument) means slot connector data is shared with the first subsolver
  DihuContext::solverStructureVisualizer()->beginChild();

  double progress = 0;
  // loop over all instances
  for (int i = 0; i < nInstancesLocal_; i++)
  {
    // output progress
    double newProgress = (double)i/nInstancesLocal_;
    if (outputInitializeThisInstance_ && this->context_.ownRankNo() == 0)
    {
      if (int(progress*10) != int(newProgress*10))
      {
        std::cout << "\b\b\b\b" << int(newProgress*100) << "%" << std::flush;
      }
    }
    progress = newProgress;

    // get the rank subset for the current instance
    std::shared_ptr<Partition::RankSubset> rankSubset = rankSubsetsLocal_[i];

    // store the rank subset containing only the own rank for the mesh of the current instance
    this->context_.partitionManager()->setRankSubsetForNextCreatedPartitioning(rankSubset);

    // call initialize on the current instance
    LOG(DEBUG) << "instance " << i << " initialize";
    instancesLocal_[i].initialize();

    // indicate that the initialize of the child is finished
    if (i == 0)
    {
      LOG(DEBUG) << " in multiple_instances, i=1, endChild and disable solverStructureVisualizer";
      DihuContext::solverStructureVisualizer()->endChild();
      DihuContext::solverStructureVisualizer()->disable();
    }
  }

  // clear rank subset for next created partitioning
  VLOG(1) << "clear rank subset for next created partitioning";
  this->context_.partitionManager()->setRankSubsetForNextCreatedPartitioning(nullptr);

  DihuContext::solverStructureVisualizer()->enable();

  // end output of progress
  if (outputInitializeThisInstance_ && this->context_.ownRankNo() == 0)
  {
    std::cout << "\b\b\b\bdone." << std::endl;
  }
  el::Loggers::addFlag(el::LoggingFlag::NewLineForContainer);

  // initialize data object with all instances
  data_.setInstancesData(instancesLocal_);

  // initialize slot connector data
  slotConnectorData_ = std::make_shared<SlotConnectorDataType>();
  slotConnectorData_->reserve(nInstancesLocal_);
  for (int i = 0; i < nInstancesLocal_; i++)
  {
    VLOG(1) << "MultipleInstances::getSlotConnectorData";
    slotConnectorData_->push_back(instancesLocal_[i].getSlotConnectorData());

    if (VLOG_IS_ON(1))
    {
      VLOG(1) << "instance " << i << "/" << nInstancesLocal_ << " is " << (*slotConnectorData_)[i];
    }
  }

  if (this->logKey_ != "")
  {
    std::stringstream logKey;
    logKey << this->logKey_ << "_init";
    Control::PerformanceMeasurement::stop(logKey.str());
  }

// #ifdef HAVE_PAT
  // PAT_region_end(1);    // end region "initialization", id 1
// #endif
}

template<typename TimeSteppingScheme>
void MultipleInstances<TimeSteppingScheme>::
run()
{
  initialize();
 
  LOG(INFO) << "MultipleInstances: " << nInstancesComputedGlobally_ << " instance" << (nInstancesComputedGlobally_ != 1? "s" : "")
    << " to be computed in total.";

#ifdef HAVE_PAT
  PAT_record(PAT_STATE_ON);
  std::string label = "computation";
  PAT_region_begin(2, label.c_str());
  LOG(INFO) << "PAT_region_begin(" << label << ")";
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
    
    // get the rank subset for the current instance
    std::shared_ptr<Partition::RankSubset> rankSubset = rankSubsetsLocal_[i];

    // store the rank subset containing only the own rank for the mesh of the current instance
    this->context_.partitionManager()->setRankSubsetForNextCreatedPartitioning(rankSubset);

    //instancesLocal_[i].reset();
    instancesLocal_[i].run();

    // avoid that solver structure file is created in every instance
    if (i == 1)
    {
      DihuContext::solverStructureVisualizer()->disable();
    }
  }

  // clear rank subset for next created partitioning
  VLOG(1) << "clear rank subset for next created partitioning";
  this->context_.partitionManager()->setRankSubsetForNextCreatedPartitioning(nullptr);

  DihuContext::solverStructureVisualizer()->enable();
  
#ifdef HAVE_PAT
  PAT_region_end(2);    // end region "computation", id 
  PAT_record(PAT_STATE_OFF);
#endif

  LOG(DEBUG) << "multipleInstances::run() complete, now call writeOutput, hasOutputWriters: "
    << std::boolalpha << this->outputWriterManager_.hasOutputWriters();

  assert(nInstancesLocal_ == instancesLocal_.size());

  // call the output writer
  if (nInstancesLocal_ > 0)
  {
    writeOwnOutput(instancesLocal_[0].numberTimeSteps(), instancesLocal_[0].endTime());
  }
  LOG(DEBUG) << "end of multiple_instances run";
}

//! return the data object
template<typename TimeSteppingScheme>
::Data::MultipleInstances<typename TimeSteppingScheme::FunctionSpace, TimeSteppingScheme> &MultipleInstances<TimeSteppingScheme>::
data()
{
  return data_;
}

template<typename TimeSteppingScheme>
void MultipleInstances<TimeSteppingScheme>::
reset()
{
  for (int i = 0; i < nInstancesLocal_; i++)
  {
    reset();
  }
}

template<typename TimeSteppingScheme>
std::shared_ptr<typename MultipleInstances<TimeSteppingScheme>::SlotConnectorDataType>
MultipleInstances<TimeSteppingScheme>::
getSlotConnectorData()
{
  // call getSlotConnectorData on all instances such that they can prepare themselves
  // (e.g. timestepping schemes call prepareForGetSlotConnectorData)
  for (int i = 0; i < nInstancesLocal_; i++)
  {
    instancesLocal_[i].getSlotConnectorData();
  }

  return slotConnectorData_;
}

template<typename TimeSteppingScheme>
std::vector<TimeSteppingScheme> &MultipleInstances<TimeSteppingScheme>::
instancesLocal()
{
  return instancesLocal_;
}

//! time of simulation
template<typename TimeSteppingScheme>
double MultipleInstances<TimeSteppingScheme>::
endTime()
{
  if (nInstancesLocal_ > 0)
    return instancesLocal_[0].endTime();
  return 0;
}

//! number of time steps in simulation time
template<typename TimeSteppingScheme>
int MultipleInstances<TimeSteppingScheme>::
numberTimeSteps()
{
  if (nInstancesLocal_ > 0)
    return instancesLocal_[0].numberTimeSteps();
  return -1;
}

template<typename TimeSteppingScheme>
void MultipleInstances<TimeSteppingScheme>::
writeOwnOutput(int timeStepNo, double currentTime, int callCountIncrement)
{
  if (nInstancesLocal_ > 0)
  {
    LOG(DEBUG) << "MultipleInstances::writeOutput, timeStepNo: " << timeStepNo << ", currentTime: " << currentTime;
    this->outputWriterManager_.writeOutput(this->data_, timeStepNo, currentTime, callCountIncrement);
  }
}

template<typename TimeSteppingScheme>
void MultipleInstances<TimeSteppingScheme>::
callOutputWriter(int timeStepNo, double currentTime, int callCountIncrement)
{
  // call the output writer of the MultipleInstances class
  writeOwnOutput(timeStepNo, currentTime, callCountIncrement);

  // call the output writers of the nested solvers
  for (int i = 0; i < nInstancesLocal_; i++)
  {
    instancesLocal_[i].callOutputWriter(timeStepNo, currentTime, callCountIncrement);
  }
}

template<typename TimeSteppingScheme>
std::string MultipleInstances<TimeSteppingScheme>::
getString(std::shared_ptr<typename MultipleInstances<TimeSteppingScheme>::SlotConnectorDataType> data)
{
  std::stringstream s;
  s << "<MultipleInstances(" << nInstancesLocal_ << "):";
  for (int i = 0; i < std::min((int)data->size(), nInstancesLocal_); i++)
  {
    if (i != 0)
      s << ", ";
    s << instancesLocal_[i].getString((*data)[i]);
  }
  s << ">";
  return s.str();
}

}  // namespace
