#include "model_order_reduction/time_stepping_scheme_ode_reduced.h"

#include <Python.h>
#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include<petscmat.h>

#include "time_stepping_scheme/time_stepping_scheme.h"

namespace ModelOrderReduction
{
template<typename TimeSteppingType>
TimeSteppingSchemeOdeReduced<TimeSteppingType>::
TimeSteppingSchemeOdeReduced(DihuContext context):
  MORBase<typename TimeSteppingType::FunctionSpace>(context["ModelOrderReduction"]),
  TimeSteppingScheme(context["ModelOrderReduction"]),
  timestepping_(context["ModelOrderReduction"]), initialized_(false)
{  
  this->specificSettings_ = this->context_.getPythonConfig();

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_, timestepping_.specificSettings());
  
  if (VLOG_IS_ON(1))
  {
    VLOG(1) << "specificSettings in TimeSteppingSchemeOdeReduced:";
    PythonUtility::printDict(this->specificSettings_);
  }
}

template<typename TimesteppingType>
std::shared_ptr<FieldVariable::FieldVariable<::FunctionSpace::Generic,1>> &TimeSteppingSchemeOdeReduced<TimesteppingType>::
solution()
{  
  return this->data_->redSolution();
}

template<typename TimeSteppingType>
void TimeSteppingSchemeOdeReduced<TimeSteppingType>::setInitialValues()
{
  PetscErrorCode ierr;
  
  Vec &solution = this->timestepping_.data().solution()->getValuesContiguous();
  Vec &redSolution=this->data_->redSolution()->getValuesContiguous();
  Mat &basisTransp = this->data_->basisTransp()->valuesGlobal();
  
  // reduction step
  ierr=MatMult(basisTransp, solution, redSolution); CHKERRV(ierr);   
}

template<typename TimeSteppingType>
void TimeSteppingSchemeOdeReduced<TimeSteppingType>::
initialize()
{
  if (initialized_)
    return;
  
  LOG(TRACE) << "TimeSteppingSchemeOdeReduced::initialize()";
  
  this->timestepping_.initialize();
  
  this->startTime_=timestepping_.startTime();
  this->endTime_=timestepping_.endTime();
  this->timeStepWidth_=timestepping_.timeStepWidth();
  this->numberTimeSteps_=timestepping_.numberTimeSteps();
  this->timeStepOutputInterval_=timestepping_.timeStepOutputInterval();
  
  LOG(DEBUG) << "timestepping_.timeStepOutputInterval() in TimeSteppingSchemeOdeReduced::initialize", timestepping_.timeStepOutputInterval();
  
  this->nReducedBases_ = PythonUtility::getOptionInt(specificSettings_, "nReducedBases", 10, PythonUtility::Positive);
  
  std::array<element_no_t, 1> nElements({this -> nReducedBases_ - 1});
  std::array<double, 1> physicalExtent({0.0});
  
  // create the functionspace for the reduced order
  std::shared_ptr<::FunctionSpace::Generic> functionSpaceRed 
    = context_.meshManager()->createFunctionSpace<::FunctionSpace::Generic>("functionSpaceReduced", nElements, physicalExtent);
  
  //assert(this->data_);
  this->data_->setFunctionSpace(functionSpaceRed);
  this->data_->setFullFunctionSpace(timestepping_.discretizableInTime().data().functionSpace());
  this->data_->initialize();
  
  MORBase<typename TimeSteppingType::FunctionSpace>::initialize();
  
  setInitialValues();
    
  initialized_=true;
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
