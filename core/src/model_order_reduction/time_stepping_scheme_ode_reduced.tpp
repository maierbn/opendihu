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
  PetscErrorCode ierr;
   
  Vec &solution = this->fullTimestepping_.data().solution()->getValuesContiguous();
  Vec &redSolution= this->data_->solution()->getValuesContiguous();
  Mat &basisTransp = this->dataMOR_->basisTransp()->valuesGlobal();
  
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
  
  this->fullTimestepping_.initialize();
   
  TimeSteppingScheme::TimeSteppingSchemeOdeBase<typename TimeSteppingType::DiscretizableInTime_Type>
  ::initialize();  
  
  if (this->specificSettingsMOR_.hasKey("nReducedBases"))
  {
    this->nReducedBases_ = this->specificSettingsMOR_.getOptionInt("nReducedBases", 10, PythonUtility::Positive);
    LOG(TRACE) << "nReducedBases: " << this->nReducedBases_;
  }
  
  if (this->specificSettingsMOR_.hasKey("nFullBases"))
  {
    this->nFullBases_ = this->specificSettingsMOR_.getOptionInt("nFullBases", 10, PythonUtility::Positive);
    LOG(TRACE) << "nFullBases: " << this->nFullBases_;
  }
  
  std::array<element_no_t, 1> nElementsRed({(this -> nReducedBases_+1)*this->discretizableInTime_.nComponents()-1});
  std::array<element_no_t, 1> nElementsFull({(this -> nFullBases_+1)*this->discretizableInTime_.nComponents()-1});
  std::array<double, 1> physicalExtent({0.0});
  
  typedef FunctionSpace::Generic GenericFunctionSpace;
  //std::shared_ptr<GenericFunctionSpace> functionSpaceRed;
  //std::shared_ptr<GenericFunctionSpace> functionSpaceFull;
  
  if(this->functionSpaceRed!=NULL)
  {
    LOG(TRACE) << "functionSpaceRed set to NULL============================";
    
    this->functionSpaceRed= NULL;
    this->functionSpaceRed= this->context_.meshManager()->template createFunctionSpace<GenericFunctionSpace>("functionSpaceReduced", nElementsRed, physicalExtent);      
  }
  else
  {
    // create the functionspace for the reduced order
    this->functionSpaceRed= this->context_.meshManager()->template createFunctionSpace<GenericFunctionSpace>("functionSpaceReduced", nElementsRed, physicalExtent);   
  
    LOG(TRACE) << "functionSpaceRed============================";
  }
  
  if(this->functionSpaceFull!=NULL)
  {
    LOG(TRACE) << "functionSpaceFull set to NULL============================";
    
    this->functionSpaceFull=nullptr;
    this->functionSpaceFull = this->context_.meshManager()->template createFunctionSpace<GenericFunctionSpace>("functionSpaceFull", nElementsFull, physicalExtent);     
  }
  else
  {
    this->functionSpaceFull = this->context_.meshManager()->template createFunctionSpace<GenericFunctionSpace>("functionSpaceFull", nElementsFull, physicalExtent);  
    
    LOG(TRACE) << "functionSpaceFull============================";
  }  
    
  this->dataMOR_->setFunctionSpace(this->functionSpaceRed);
  this->dataMOR_->setFullFunctionSpace(this->functionSpaceFull);
  //this->dataMOR_->setFullFunctionSpace(this->fullTimestepping_.discretizableInTime_.functionSpace());
  
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
