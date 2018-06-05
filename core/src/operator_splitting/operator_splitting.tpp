#include "operator_splitting/operator_splitting.h"

#include "utility/python_utility.h"
#include "data_management/time_stepping.h"

namespace OperatorSplitting
{

template<typename TimeStepping1, typename TimeStepping2>
OperatorSplitting<TimeStepping1, TimeStepping2>::
OperatorSplitting(DihuContext context, std::string schemeName) :
  ::TimeSteppingScheme::TimeSteppingScheme(context),
  timeStepping1_(context_[schemeName]["Term1"]),
  timeStepping2_(context_[schemeName]["Term2"]), initialized_(false)
{
  PyObject *topLevelSettings = context_.getPythonConfig();
  specificSettings_ = PythonUtility::getOptionPyObject(topLevelSettings, schemeName);
  outputWriterManager_.initialize(specificSettings_);
}

template<typename TimeStepping1, typename TimeStepping2>
void OperatorSplitting<TimeStepping1, TimeStepping2>::
initialize()
{
  if (initialized_)
    return;
  LOG(TRACE) << "  OperatorSplitting::initialize";

  TimeSteppingScheme::initialize();

  LOG(TRACE) << "  OperatorSplitting::initialize done, timeSpan=[" <<this->startTime_<<","<<this->endTime_<<"]"
    <<", n steps: "<<this->numberTimeSteps_;

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

  outputData1_ = PythonUtility::getOptionBool(specificSettings_, "outputData1", true);
  outputData2_ = PythonUtility::getOptionBool(specificSettings_, "outputData2", true);
  
  initialized_ = true;
}

template<typename TimeStepping1, typename TimeStepping2>
Vec &OperatorSplitting<TimeStepping1, TimeStepping2>::
solution()
{
  return timeStepping1_.solution();
}

template<typename TimeStepping1, typename TimeStepping2>
void OperatorSplitting<TimeStepping1, TimeStepping2>::
run()
{
  // initialize data structurures
  initialize();

  // run simulation
  advanceTimeSpan();
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


};    // namespace