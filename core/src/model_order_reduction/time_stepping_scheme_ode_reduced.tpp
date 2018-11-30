#include "model_order_reduction/time_stepping_scheme_ode_reduced.h"

#include <Python.h>
#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include<petscmat.h>
#include "mesh/mesh_manager.h"
#include "function_space/function_space.h"

namespace ModelOrderReduction
{
  template<typename TimeSteppingType>
TimeSteppingSchemeOdeReduced<TimeSteppingType>::
TimeSteppingSchemeOdeReduced(DihuContext context, std::string name):
  MORBase<typename TimeSteppingType::FunctionSpace>(context["ModelOrderReduction"]),
  TimeSteppingScheme::TimeSteppingSchemeOdeBase<typename TimeSteppingType::DiscretizableInTime_Type>(context["ModelOrderReduction"],name),
  timestepping_(context["ModelOrderReduction"]), initialized_(false)
{ /* 
  this->specificSettings_ = this->context_.getPythonConfig();

  // initialize output writers
  this->outputWriterManager_.initialize(this->context_, timestepping_.specificSettings());
  
  if (VLOG_IS_ON(1))
  {
    VLOG(1) << "specificSettings in TimeSteppingSchemeOdeReduced:";
    PythonUtility::printDict(this->specificSettings_.pyObject());
  }
  */
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
  /*
  PetscErrorCode ierr;
   
  Vec &solution = this->discretizableInTime_.data().solution()->getValuesContiguous();
  Vec &redSolution=this->data_->solution()->getValuesContiguous();
  Mat &basisTransp = this->dataMOR_->basisTransp()->valuesGlobal();
  
  // reduction step
  ierr=MatMult(basisTransp, solution, redSolution); CHKERRV(ierr); 
  */
}

template<typename TimeSteppingType>
void TimeSteppingSchemeOdeReduced<TimeSteppingType>::
initialize()
{
  if (initialized_)
    return;
  
  TimeSteppingScheme::TimeSteppingSchemeOdeBase<typename TimeSteppingType::DiscretizableInTime_Type>::initialize();
  
  LOG(TRACE) << "TimeSteppingSchemeOdeReduced::initialize()";
  
  this->timestepping_.initialize();
  
  this->startTime_=timestepping_.startTime();
  this->endTime_=timestepping_.endTime();
  this->timeStepWidth_=timestepping_.timeStepWidth();
  this->numberTimeSteps_=timestepping_.numberTimeSteps();
  this->timeStepOutputInterval_=timestepping_.timeStepOutputInterval();
  
  LOG(DEBUG) << "timestepping_.timeStepOutputInterval() in TimeSteppingSchemeOdeReduced::initialize", timestepping_.timeStepOutputInterval();
  
  this->nReducedBases_ = this->specificSettings_.getOptionInt("nReducedBases", 10, PythonUtility::Positive);
  
  std::array<element_no_t, 1> nElements({this -> nReducedBases_ - 1});
  std::array<double, 1> physicalExtent({0.0});
  
  typedef FunctionSpace::Generic GenericFunctionSpace;
  
  // create the functionspace for the reduced order
  std::shared_ptr<GenericFunctionSpace> functionSpaceRed 
    = this->context_.meshManager()->template createFunctionSpace<GenericFunctionSpace>("functionSpaceReduced", nElements, physicalExtent);
  
  //assert(this->data_);
  this->data_->setFunctionSpace(functionSpaceRed);
  this->dataMOR_->setFullFunctionSpace(this->discretizableInTime_.functionSpace());
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
bool TimeSteppingSchemeOdeReduced<TimeSteppingType>::
knowsMeshType()
{
  return false;
}
  
} //namespace
