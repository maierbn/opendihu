#include "control/multiple_instances.h"

#include "data_management/multiple_instances.h"

namespace Control
{

template<class TimeSteppingScheme>
MultipleInstances<TimeSteppingScheme>::
MultipleInstances(DihuContext context) :
  context_(context["MultipleInstances"]), data_(context_)
{
  specificSettings_ = context_.getPythonConfig();
  outputWriterManager_.initialize(specificSettings_);
  
  // extract the number of instances
  nInstances_ = PythonUtility::getOptionInt(specificSettings_, "nInstances", 1, PythonUtility::Positive);
   
  // get the config for the firsetInstancesDatast instance from the list
  PyObject *instanceConfig = PythonUtility::getOptionListBegin<PyObject *>(specificSettings_, "instances");

  // loop over other entries of list
  
  // create all instances
  for (int i = 0; i < nInstances_; i++)
  {
    if (PythonUtility::getOptionListEnd(specificSettings_, "instances"))
    {
      LOG(ERROR) << "Could only create " << i << " instances from the given instances config python list, but nInstances = " << nInstances_;
      nInstances_ = i;
      break;
    }
   
    instances_.emplace_back(context_.createSubContext(instanceConfig));
    
    // extract config for next instance
    PythonUtility::getOptionListNext<PyObject *>(specificSettings_, "instances", instanceConfig);
  }
  
  if (!PythonUtility::getOptionListEnd(specificSettings_, "instances"))
  {
    LOG(WARNING) << "Only " << nInstances_ << " instances were created, but more configurations are given.";
  }
}

template<class TimeSteppingScheme>
void MultipleInstances<TimeSteppingScheme>::
advanceTimeSpan()
{
  // This method advances the simulation by the specified time span. It will be needed when this MultipleInstances object is part of a parent control element, like a coupling to 3D model.
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
  
  data_.setInstancesData(instances_);
}

template<class TimeSteppingScheme>
void MultipleInstances<TimeSteppingScheme>::
run()
{
  initialize();
 
  #pragma omp parallel for
  for (int i = 0; i < nInstances_; i++)
  {
    instances_[i].run();
  }
  
  this->outputWriterManager_.writeOutput(this->data_);
}

template<class TimeSteppingScheme>
bool MultipleInstances<TimeSteppingScheme>::
knowsMeshType()
{
  // This is a dummy method that is currently not used, it is only important if we want to map between multiple data sets.
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

};