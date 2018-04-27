#include "control/multiple_instances.h"

namespace Control
{

template<class TimeSteppingScheme>
MultipleInstances<TimeSteppingScheme>::
MultipleInstances(DihuContext context) : 
  context_(context["MultipleInstances"])
{
  // extract the number of instances
  nInstances_ = PythonUtility::getOptionInt(context_.getPythonConfig(), "nInstances", 1, PythonUtility::Positive);
  
  // create all instances
  for (int i = 0; i < nInstances_; i++)
  {
    instances_.emplace_back(context_);
  }
}

template<class TimeSteppingScheme>
void MultipleInstances<TimeSteppingScheme>::
advanceTimeSpan()
{
  for (int i = 0; i < nInstances_; i++)
  {
    instances_[i].advanceTimeSpan();
  }
}

template<class TimeSteppingScheme>
void MultipleInstances<TimeSteppingScheme>::
setTimeSpan(double startTime, double endTime)
{
  for (int i = 0; i < nInstances_; i++)
  {
    instances_[i].setTimeSpan(startTime, endTime);
  }
}

template<class TimeSteppingScheme>
void MultipleInstances<TimeSteppingScheme>::
initialize()
{
  for (int i = 0; i < nInstances_; i++)
  {
    instances_[i].initialize();
  }
}

template<class TimeSteppingScheme>
void MultipleInstances<TimeSteppingScheme>::
run()
{
  for (int i = 0; i < nInstances_; i++)
  {
    instances_[i].run();
  }
}

template<class TimeSteppingScheme>
bool MultipleInstances<TimeSteppingScheme>::
knowsMeshType()
{
  assert(nInstances_ > 0);
  return instances_[0].knowsMeshType();
}

template<class TimeSteppingScheme>
Vec &MultipleInstances<TimeSteppingScheme>::
solution()
{
  assert(nInstances_ > 0);
  return instances_[0].solution();
}

template<class TimeSteppingScheme>
SolutionVectorMapping &MultipleInstances<TimeSteppingScheme>::
solutionVectorMapping()
{
  assert(nInstances_ > 0);
  return instances_[0].solutionVectorMapping();
}

};