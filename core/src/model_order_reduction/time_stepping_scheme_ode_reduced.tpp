#pragma once

#include "model_order_reduction/time_stepping_scheme_ode_reduced.h"

#include <Python.h>
#include "utility/python_utility.h"
#include "time_stepping_scheme/time_stepping_scheme.h"

namespace ModelOrderReduction
{
  template<typename TimeSteppingType>
  TimeSteppingSchemeOdeReduced<TimeSteppingType>::
  TimeSteppingSchemeOdeReduced(DihuContext context):
  MORBase<FunctionSpace::Generic>(context["ModelOrderReduction"]), TimeSteppingScheme(context["ModelOrderReduction"]),
  timestepping_(context_["ModelOrderReduction"]), solutionVectorMapping_(SolutionVectorMapping()), initialized_(false)
  {
    PyObject *topLevelSettings = this->context_.getPythonConfig();
    this->specificSettings_ = PythonUtility::getOptionPyObject(topLevelSettings, "ModelOrderReduction");
  }
  
  template<typename TimesteppingType>
  FieldVariable::FieldVariable<FunctionSpace::Generic,1> &TimeSteppingSchemeOdeReduced<TimesteppingType>::
  solution()
  {
    // invert the PETSc Vec object
    //std::shared_ptr<Partition::MeshPartition<FunctionSpace::Generic>> partition = this->functionSpace_->meshPartition();
    //creating the ptr each time, to improve
    //return std::make_shared<PartitionedPetscVec<FunctionSpace::Generic>>(partition, this->data_->Solution(), "redSolution");
    
    return;
  }
  
  template<typename TimeSteppingType>
  SolutionVectorMapping &TimeSteppingSchemeOdeReduced<TimeSteppingType>::
  solutionVectorMapping()
  {
    return solutionVectorMapping_;
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
  
template<typename TimeSteppingType>
void TimeSteppingSchemeOdeReduced<TimeSteppingType>::
reset()
{
  TimeSteppingScheme::reset();
  timestepping_.reset();
  
  initialized_ = false;
}
  
  template<typename TimeSteppingType>
  bool TimeSteppingSchemeOdeReduced<TimeSteppingType>::
  knowsMeshType()
  {
    return false;
  }
  
} //namespace