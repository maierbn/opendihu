#include "operator_splitting/operator_splitting.h"

#ifdef HAVE_PAT
#include <pat_api.h>    // perftools, only available on hazel hen
#endif

#include "utility/python_utility.h"
#include "data_management/time_stepping/time_stepping.h"
#include "control/diagnostic_tool/performance_measurement.h"
#include "control/diagnostic_tool/solver_structure_visualizer.h"
#include "mesh/mesh_manager/mesh_manager.h"

namespace OperatorSplitting
{

template<typename TimeStepping1, typename TimeStepping2>
OperatorSplitting<TimeStepping1, TimeStepping2>::
OperatorSplitting(DihuContext context, std::string schemeName) :
  ::TimeSteppingScheme::TimeSteppingScheme(context),
  timeStepping1_(context_[schemeName]["Term1"]),
  timeStepping2_(context_[schemeName]["Term2"]),
  initialized_(false)
{

  PythonConfig topLevelSettings = context_.getPythonConfig();
  specificSettings_ = PythonConfig(topLevelSettings, schemeName);
  schemeName_ = schemeName;

  outputConnection_ = std::make_shared<OutputConnection>(specificSettings_);
}

template<typename TimeStepping1, typename TimeStepping2>
void OperatorSplitting<TimeStepping1, TimeStepping2>::
reset()
{
  initialized_ = false;
  timeStepping1_.reset();
  timeStepping2_.reset();
}

template<typename TimeStepping1, typename TimeStepping2>
void OperatorSplitting<TimeStepping1, TimeStepping2>::
initialize()
{
  if (initialized_)
    return;
  LOG(TRACE) << "  OperatorSplitting::initialize";

  // add this solver to the solvers diagram, which is a SVG file that will be created at the end of the simulation.
  DihuContext::solverStructureVisualizer()->addSolver(schemeName_, true);   // hasInternalConnectionToFirstNestedSolver=true (the last argument) means output connector data is shared with the first subsolver

  TimeSteppingScheme::initialize();
  timeStepOutputInterval_ = specificSettings_.getOptionInt("timeStepOutputInterval", 100, PythonUtility::Positive);

  if (specificSettings_.hasKey("transferSlotName"))
    LOG(WARNING) << "In " << specificSettings_ << ", option \"transferSlotName\" is no longer used, use \"connectedSlotsTerm1To2\" and \"connectedSlotsTerm2To1\" instead.";

  LOG(TRACE) << "  OperatorSplitting::initialize done, timeSpan=[" << this->startTime_<< "," << this->endTime_<< "]"
    << ", n steps: " << this->numberTimeSteps_;


  // indicate in solverStructureVisualizer that now a child solver will be initialized
  DihuContext::solverStructureVisualizer()->beginChild("Term1");

  // initialize time stepping objects
  LOG(DEBUG) << "  OperatorSplitting::initialize timeStepping1";
  timeStepping1_.initialize();

  // indicate in solverStructureVisualizer that the child solver initialization is done
  DihuContext::solverStructureVisualizer()->endChild();

  // indicate in solverStructureVisualizer that now a child solver will be initialized
  DihuContext::solverStructureVisualizer()->beginChild("Term2");

  LOG(DEBUG) << "  OperatorSplitting::initialize timeStepping2";
  timeStepping2_.initialize();

  // indicate in solverStructureVisualizer that the child solver initialization is done
  DihuContext::solverStructureVisualizer()->endChild();

  LOG(DEBUG) << "initialize mappings between meshes \"" << timeStepping1_.data().functionSpace()->meshName() << "\" and \""
    << timeStepping2_.data().functionSpace()->meshName() << "\".";
  context_.meshManager()->initializeMappingsBetweenMeshes<typename TimeStepping1::FunctionSpace,typename TimeStepping2::FunctionSpace>(
    timeStepping1_.data().functionSpace(), timeStepping2_.data().functionSpace());

  // log endTime parameters
  Control::PerformanceMeasurement::setParameter("endTime", endTime_);

  // compose logging keys
  logKeyTimeStepping1AdvanceTimeSpan_ = this->durationLogKey_ + std::string("_advanceTimeSpan1");  ///< key for logging of the duration of the advanceTimeSpan() call of timeStepping1
  logKeyTimeStepping2AdvanceTimeSpan_ = this->durationLogKey_ + std::string("_advanceTimeSpan2");  ///< key for logging of the duration of the advanceTimeSpan() call of timeStepping2
  logKeyTransfer12_ = this->durationLogKey_ + std::string("_transfer12");  ///< key for logging of the duration of data transfer from timestepping 1 to 2
  logKeyTransfer21_ = this->durationLogKey_ + std::string("_transfer21");  ///< key for logging of the duration of data transfer from timestepping 2 to 1

  // set the outputConnection information about how the slots are connected to appear in the solver diagram
  DihuContext::solverStructureVisualizer()->addOutputConnection(outputConnection_);

  // set the outputConnectorData for the solverStructureVisualizer to appear in the solver diagram
  DihuContext::solverStructureVisualizer()->setOutputConnectorData(getOutputConnectorData());

  initialized_ = true;
}

template<typename TimeStepping1, typename TimeStepping2>
void OperatorSplitting<TimeStepping1, TimeStepping2>::
setRankSubset(Partition::RankSubset rankSubset)
{
  timeStepping1_.setRankSubset(rankSubset);
  timeStepping2_.setRankSubset(rankSubset);
}

  
template<typename TimeStepping1, typename TimeStepping2>
void OperatorSplitting<TimeStepping1, TimeStepping2>::
run()
{
  // initialize data structurures
  initialize();

#ifdef HAVE_PAT
  PAT_record(PAT_STATE_ON);
  std::string label = "computation";
  PAT_region_begin(2, label.c_str());
  LOG(INFO) << "PAT_region_begin(" << label << ")";
#endif

  // run simulation
  advanceTimeSpan();

#ifdef HAVE_PAT
  PAT_region_end(2);    // end region "computation", id 
  PAT_record(PAT_STATE_OFF);
#endif
}

template<typename TimeStepping1, typename TimeStepping2>
std::shared_ptr<typename OperatorSplitting<TimeStepping1, TimeStepping2>::OutputConnectorDataType>
OperatorSplitting<TimeStepping1, TimeStepping2>::
getOutputConnectorData()
{
  return timeStepping1_.getOutputConnectorData();
}

template<typename TimeStepping1, typename TimeStepping2>
typename OperatorSplitting<TimeStepping1, TimeStepping2>::Data &OperatorSplitting<TimeStepping1, TimeStepping2>::
data()
{
  return timeStepping1_.data();
}

//! get a reference to the first timestepping object
template<typename TimeStepping1, typename TimeStepping2>
TimeStepping1 &OperatorSplitting<TimeStepping1, TimeStepping2>::
timeStepping1()
{
  return timeStepping1_;
}

//! get a reference to the second timestepping object
template<typename TimeStepping1, typename TimeStepping2>
TimeStepping2 &OperatorSplitting<TimeStepping1, TimeStepping2>::
timeStepping2()
{
  return timeStepping2_;
}

//! output the given data for debugging
template<typename TimeStepping1, typename TimeStepping2>
std::string OperatorSplitting<TimeStepping1, TimeStepping2>::
getString(std::shared_ptr<typename OperatorSplitting<TimeStepping1, TimeStepping2>::OutputConnectorDataType> data)
{
  std::stringstream s;
  s << "<" << schemeName_ << ",Term1:" << data << ">";
  return s.str();
}


}  // namespace
