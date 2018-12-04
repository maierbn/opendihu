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
  fullTimestepping_(context["ModelOrderReduction"]), initialized_(false)
{  
  this->specificSettingsMOR_ = this->context_.getPythonConfig(); 
  this->data_ = std::make_shared <Data::TimeStepping<typename TimeSteppingType::FunctionSpace, TimeSteppingType::DiscretizableInTime_Type::nComponents()>>(context); // create data object
  
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
  
  LOG(TRACE) << "TimeSteppingSchemeOdeReduced::initialize()";
  
  this->fullTimestepping_.initialize();
  
  LOG(TRACE) << "Finished fullTimestepping_.initialize()==============================";
  
  
  TimeSteppingScheme::TimeSteppingSchemeOdeBase<typename TimeSteppingType::DiscretizableInTime_Type>
  ::initialize();  
  
  LOG(TRACE) << "Finished TimeSteppingSchemeOdeBase.initialize()==============================";
  
  if (this->specificSettingsMOR_.hasKey("nReducedBases"))
  {
    this->nReducedBases_ = this->specificSettingsMOR_.getOptionInt("nReducedBases", 10, PythonUtility::Positive);
    LOG(TRACE) << "nReducedBases: " << this->nReducedBases_;
  }
  
  std::array<element_no_t, 1> nElements({(this -> nReducedBases_ - 1)*this->discretizableInTime_.nComponents()});
  std::array<double, 1> physicalExtent({0.0});
  
  typedef FunctionSpace::Generic GenericFunctionSpace;
  
  // create the functionspace for the reduced order
  std::shared_ptr<GenericFunctionSpace> functionSpaceRed 
    = this->context_.meshManager()->template createFunctionSpace<GenericFunctionSpace>("functionSpaceReduced", nElements, physicalExtent);  
   
  //TimeSteppingScheme::TimeSteppingSchemeOdeBase<typename TimeSteppingType::DiscretizableInTime_Type>::initialize(functionSpaceRed);  
  //assert(functionSpaceRed->meshPartition());
  //assert(this->data_);
  //assert(TimeSteppingScheme::TimeSteppingSchemeOdeBase<typename TimeSteppingType::DiscretizableInTime_Type>::data_);
  //this->data_->setFunctionSpace(functionSpaceRed);
  
  //this->data_->initialize();
    
  this->dataMOR_->setFunctionSpace(functionSpaceRed);
  this->dataMOR_->setFullFunctionSpace(this->discretizableInTime_.functionSpace());
  MORBase<typename TimeSteppingType::FunctionSpace>::initialize();  
  
  //setInitialValues();
    
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
