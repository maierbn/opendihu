#include "model_order_reduction/time_stepping_scheme_ode_reduced.h"

#include <Python.h>
#include "utility/python_utility.h"
#include "utility/petsc_utility.h"
#include<petscmat.h>
#include "mesh/mesh_manager.h"
#include "function_space/function_space.h"
#include "time_stepping_scheme/time_stepping_scheme.h"
#include "time_stepping_scheme/time_stepping_scheme_ode.h"
#include "data_management/time_stepping/time_stepping.h"

namespace ModelOrderReduction
{
  template<typename TimeSteppingType>
TimeSteppingSchemeOdeReduced<TimeSteppingType>::
TimeSteppingSchemeOdeReduced(DihuContext context, std::string name):
  MORBase<typename TimeSteppingType::FunctionSpace>(context["ModelOrderReduction"]),
  ::TimeSteppingScheme::TimeSteppingSchemeOdeBase<::FunctionSpace::Generic,1>(context["ModelOrderReduction"],name),
  fullTimestepping_(context["ModelOrderReduction"]), initialized_(false)
{  
  this->specificSettingsMOR_ = this->context_.getPythonConfig();
  
  if (this->specificSettingsMOR_.hasKey("nReducedBases"))
  {
    this->nReducedBases_ = this->specificSettingsMOR_.getOptionInt("nReducedBases", 10, PythonUtility::Positive);
    LOG(TRACE) << "nReducedBases: " << this->nReducedBases_;
  }
  
  if (this->specificSettingsMOR_.hasKey("nRowsSnapshots"))
  {
    this->nRowsSnapshots_ = this->specificSettingsMOR_.getOptionInt("nRowsSnapshots", 10, PythonUtility::Positive);
    LOG(TRACE) << "nRowsSnapshots: " << this->nRowsSnapshots_;
  }
  
  std::array<element_no_t, 1> nElementsRed({this -> nReducedBases_});
  std::array<element_no_t, 1> nElementsRows({this -> nRowsSnapshots_});
  std::array<double, 1> physicalExtent({0.0}); 
  
  typedef ::FunctionSpace::Generic GenericFunctionSpace;
  
  if(this->context_.meshManager()->hasFunctionSpace("functionSpaceReduced"))
  {
    LOG(TRACE) << "functionSpaceRed set to NULL============================";
    this->functionSpaceRed= this->context_.meshManager()->template functionSpace<GenericFunctionSpace>("functionSpaceReduced");
    //this->functionSpaceRed= NULL;
    //this->functionSpaceRed= this->context_.meshManager()->template createFunctionSpace<GenericFunctionSpace>("functionSpaceReduced", nElementsRed, physicalExtent);      
  }
  else
  {
    // create the functionspace for the reduced order
    this->functionSpaceRed= this->context_.meshManager()->template createFunctionSpace<GenericFunctionSpace>("functionSpaceReduced", nElementsRed, physicalExtent);   
    
    LOG(TRACE) << "functionSpaceRed============================";
  }
  
  if(this->context_.meshManager()->hasFunctionSpace("functionSpaceRowsSnapshots"))
  {
    LOG(TRACE) << "functionSpaceRowsSnapshots set to NULL============================";
    this->functionSpaceRowsSnapshots= this->context_.meshManager()->template functionSpace<GenericFunctionSpace>("functionSpaceRowsSnapshots");
    
    //this->functionSpaceRowsSnapshots=nullptr;
    //this->functionSpaceRowsSnapshots = this->context_.meshManager()->template createFunctionSpace<GenericFunctionSpace>("functionSpaceFull", nElementsFull, physicalExtent);     
  }
  else
  {
    this->functionSpaceRowsSnapshots = this->context_.meshManager()->template createFunctionSpace<GenericFunctionSpace>("functionSpaceRowsSnapshots", nElementsRows, physicalExtent);  
    
    LOG(TRACE) << "functionSpaceRowsSnapshots============================";
  }
  
  this->data_ = std::make_shared <::Data::TimeStepping<::FunctionSpace::Generic,1>>(context); // create data object

}

//template<typename TimesteppingType>
//std::shared_ptr<FieldVariable::FieldVariable<::FunctionSpace::Generic,1>> &TimeSteppingSchemeOdeReduced<TimesteppingType>::
//solution()
//{  
//  return this->data_->solution();
//}

template<typename TimeSteppingType>
void TimeSteppingSchemeOdeReduced<TimeSteppingType>::setInitialValues()
{  
  PetscErrorCode ierr;
   
  Vec &solution = this->fullTimestepping_.data().solution()->getValuesContiguous();
  Vec &redSolution= this->data().solution()->getValuesContiguous();
  Mat &basisTransp = this->dataMOR_->basisTransp()->valuesGlobal();
  
  PetscInt mat_sz_1, mat_sz_2;
  PetscInt solution_size, redSolution_size;
  
  VecGetSize(solution,&solution_size);
  VecGetSize(redSolution,&redSolution_size);
  MatGetSize(basisTransp,&mat_sz_1,&mat_sz_2);
  
  LOG(DEBUG) << "setInitialValues() solution_size: " << solution_size << "========================";
  LOG(DEBUG) << "setInitialValues() redSolution_size: " << redSolution_size; 
  LOG(DEBUG) << "setInitialValues() mat_sz_1: " << mat_sz_1 << " mat_sz_2: " << mat_sz_2 << "========================";;
  
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
   
  ::TimeSteppingScheme::TimeSteppingSchemeOdeBase<::FunctionSpace::Generic,1>::initialize(); 
  
  //this->data_->setOutputComponentNo(outputComponentNo);
  
  /*
  if (this->specificSettingsMOR_.hasKey("nReducedBases"))
  {
    this->nReducedBases_ = this->specificSettingsMOR_.getOptionInt("nReducedBases", 10, PythonUtility::Positive);
    LOG(TRACE) << "nReducedBases: " << this->nReducedBases_;
  }
 
  if (this->specificSettingsMOR_.hasKey("nRowsSnapshots"))
  {
    this->nFullBases_ = this->specificSettingsMOR_.getOptionInt("nRowsSnapshots", 10, PythonUtility::Positive);
    LOG(TRACE) << "nRowsSnapshots: " << this->nRowsSnapshots_;
  }
  
  std::array<element_no_t, 1> nElementsRed({this -> nReducedBases_});
  std::array<element_no_t, 1> nElementsRows({this -> nRowsSnapshots_});
  std::array<double, 1> physicalExtent({0.0});
  
  typedef FunctionSpace::Generic GenericFunctionSpace;

  if(this->context_.meshManager()->hasFunctionSpace("functionSpaceReduced"))
  {
    LOG(TRACE) << "functionSpaceRed set to NULL============================";
    this->functionSpaceRed= this->context_.meshManager()->template functionSpace<GenericFunctionSpace>("functionSpaceReduced");
    //this->functionSpaceRed= NULL;
    //this->functionSpaceRed= this->context_.meshManager()->template createFunctionSpace<GenericFunctionSpace>("functionSpaceReduced", nElementsRed, physicalExtent);      
  }
  else
  {
    // create the functionspace for the reduced order
    this->functionSpaceRed= this->context_.meshManager()->template createFunctionSpace<GenericFunctionSpace>("functionSpaceReduced", nElementsRed, physicalExtent);   
  
    LOG(TRACE) << "functionSpaceRed============================";
  }
  
  if(this->context_.meshManager()->hasFunctionSpace("functionSpaceRowsSnapshots"))
  {
    LOG(TRACE) << "functionSpaceRowsSnapshots set to NULL============================";
    this->functionSpaceRowsSnapshots= this->context_.meshManager()->template functionSpace<GenericFunctionSpace>("functionSpaceRowsSnapshots");
    
    //this->functionSpaceRowsSnapshots=nullptr;
    //this->functionSpaceRowsSnapshots = this->context_.meshManager()->template createFunctionSpace<GenericFunctionSpace>("functionSpaceFull", nElementsFull, physicalExtent);     
  }
  else
  {
    this->functionSpaceRowsSnapshots = this->context_.meshManager()->template createFunctionSpace<GenericFunctionSpace>("functionSpaceRowsSnapshots", nElementsFull, physicalExtent);  
    
    LOG(TRACE) << "functionSpaceRowsSnapshots============================";
  }  
  */

  this->dataMOR_->setFunctionSpace(this->functionSpaceRed);
  this->dataMOR_->setFunctionSpaceRows(this->functionSpaceRowsSnapshots);
  
  assert(functionSpaceRowsSnapshots->meshPartition());   // assert that the function space was retrieved correctly
  data_->setFunctionSpace(functionSpaceRowsSnapshots);
  this->data().setOutputComponentNo(0);
  this->data_->initialize();
  
  MORBase<typename TimeSteppingType::FunctionSpace>::initialize();  
  
  setInitialValues(); //necessary for the implicit scheme
    
  initialized_=true;
}
/*
template<int nStates, typename FunctionSpaceType>
void TimeSteppingSchemeOdeReduced<CellmlAdapter<nStates, FunctionSpaceType>>::
initialize()
{
  
}
*/
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
