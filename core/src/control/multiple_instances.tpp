#include "control/multiple_instances.h"

#include <omp.h>

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
  
  //VLOG(1) << "MultipleInstances constructor, settings: " << specificSettings_;
  
  // extract the number of instances
  nInstances_ = PythonUtility::getOptionInt(specificSettings_, "nInstances", 1, PythonUtility::Positive);
   
  // parse all instance configs 
  std::vector<PyObject *> instanceConfigs;
  
  // get the config for the firsetInstancesDatast instance from the list
  PyObject *instanceConfig = PythonUtility::getOptionListBegin<PyObject *>(specificSettings_, "instances");

  int i = 0;
  for(;
      !PythonUtility::getOptionListEnd(specificSettings_, "instances") && i < nInstances_; 
      PythonUtility::getOptionListNext<PyObject *>(specificSettings_, "instances", instanceConfig), i++)
  {
    instanceConfigs.push_back(instanceConfig);
    VLOG(1) << "i = " << i << ", instanceConfig = " << instanceConfig;
  }
    
  if (i < nInstances_)
  {
    LOG(ERROR) << "Could only create " << i << " instances from the given instances config python list, but nInstances = " << nInstances_;
  }
    
  if (!PythonUtility::getOptionListEnd(specificSettings_, "instances"))
  {
    LOG(ERROR) << "Only " << nInstances_ << " instances were created, but more configurations are given.";
  }
   
  // create all instances
  for (PyObject *instanceConfig: instanceConfigs)
  {
    VLOG(1) << "create sub context";
    instances_.emplace_back(context_.createSubContext(instanceConfig));
    
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
    if (omp_get_thread_num() == 0)
      LOG(INFO) << "running with " << omp_get_num_threads() << " OpenMP threads";
    
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