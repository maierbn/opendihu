#include "model_order_reduction/time_stepping_scheme_ode_reduced.h"

#include <Python.h>
#include "utility/python_utility.h"
#include "time_stepping_scheme/time_stepping_scheme.h"

namespace ModelOrderReduction
{
  template<typename TimeSteppingType>
  TimeSteppingSchemeOdeReduced<TimeSteppingType>::
  TimeSteppingSchemeOdeReduced(DihuContext context):
  MORBase(context["ModelOrderReduction"]), TimeSteppingScheme(context["ModelOrderReduction"]),
  timestepping_(context_["ModelOrderReduction"]), initialized_(false) 
  {
    PyObject *topLevelSettings = this->context_.getPythonConfig();
    this->specificSettings_ = PythonUtility::getOptionPyObject(topLevelSettings, "ModelOrderReduction");
  }
  
  template<typename TimeSteppingType>
  void TimeSteppingSchemeOdeReduced<TimeSteppingType>::
  initialize()
  {
    if (initialized_)
      return;
    
    LOG(TRACE) << "TimeSteppingSchemeOdeReduced::initialize()";
    
    MORBase::initialize();
    this->timestepping_.initialize();   
  }
  
  template<typename TimeSteppingType>
  void TimeSteppingSchemeOdeReduced<TimeSteppingType>::
  run()
  {
    // initialize
    this->initialize();
    
    // do simulations
    this->advanceTimeSpan();    
  }
  
} //namespace