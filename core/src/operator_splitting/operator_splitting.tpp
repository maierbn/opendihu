#include "operator_splitting/operator_splitting.h"

#ifdef HAVE_PAT
#include <pat_api.h>    // perftools, only available on hazel hen
#endif

#include "utility/python_utility.h"
#include "data_management/time_stepping/time_stepping.h"
#include "control/performance_measurement.h"
#include "mesh/mesh_manager.h"

namespace OperatorSplitting
{

template<typename TimeStepping1, typename TimeStepping2>
OperatorSplitting<TimeStepping1, TimeStepping2>::
OperatorSplitting(DihuContext context, std::string schemeName) :
  ::TimeSteppingScheme::TimeSteppingScheme(context),
  timeStepping1_(context_[schemeName]["Term1"]),
  timeStepping2_(context_[schemeName]["Term2"]), initialized_(false)
{

  PythonConfig topLevelSettings = context_.getPythonConfig();
  specificSettings_ = PythonConfig(topLevelSettings, schemeName);
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

  TimeSteppingScheme::initialize();
  timeStepOutputInterval_ = specificSettings_.getOptionInt("timeStepOutputInterval", 100, PythonUtility::Positive);

  LOG(TRACE) << "  OperatorSplitting::initialize done, timeSpan=[" << this->startTime_<< "," << this->endTime_<< "]"
    << ", n steps: " << this->numberTimeSteps_;

  // initialize time stepping objects, if only one knows its MeshType, initialize that first
  //(e.g. CellML-adapter does not know, because it is independent of the type of a mesh)
  if (timeStepping2_.knowsMeshType() && !timeStepping1_.knowsMeshType())
  {
    LOG(DEBUG) << "  OperatorSplitting::initialize timeStepping2";
    timeStepping2_.initialize();
    LOG(DEBUG) << "  OperatorSplitting::initialize timeStepping1";
    timeStepping1_.initialize();
  }
  else
  {
    LOG(DEBUG) << "  OperatorSplitting::initialize timeStepping1";
    timeStepping1_.initialize();
    LOG(DEBUG) << "  OperatorSplitting::initialize timeStepping2";
    timeStepping2_.initialize();
  }

  context_.meshManager()->initializeMappingsBetweenMeshes<typename TimeStepping1::FunctionSpace,typename TimeStepping2::FunctionSpace>(
    timeStepping1_.data().functionSpace(), timeStepping2_.data().functionSpace());

  // log endTime parameters
  Control::PerformanceMeasurement::setParameter("endTime", endTime_);

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
typename OperatorSplitting<TimeStepping1, TimeStepping2>::TransferableSolutionDataType OperatorSplitting<TimeStepping1, TimeStepping2>::
getSolutionForTransfer()
{
  return timeStepping2_.getSolutionForTransfer();
}

template<typename TimeStepping1, typename TimeStepping2>
bool OperatorSplitting<TimeStepping1, TimeStepping2>::
knowsMeshType()
{
  return timeStepping1_.knowsMeshType() && timeStepping2_.knowsMeshType();
}

template<typename TimeStepping1, typename TimeStepping2>
typename OperatorSplitting<TimeStepping1, TimeStepping2>::Data &OperatorSplitting<TimeStepping1, TimeStepping2>::
data()
{
  return timeStepping1_.data();
}


}  // namespace
