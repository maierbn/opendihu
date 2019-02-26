#include "model_order_reduction/time_stepping_reduced_explicit.h"

#include <Python.h>
#include<petscmat.h>
#include "utility/python_utility.h"
#include "utility/petsc_utility.h"


namespace ModelOrderReduction
{
  template<typename TimeSteppingExplicitType>
  TimeSteppingSchemeOdeReducedExplicit<TimeSteppingExplicitType>::
  TimeSteppingSchemeOdeReducedExplicit(DihuContext context,std::string name):
  TimeSteppingSchemeOdeReduced<TimeSteppingExplicitType>(context,name)
  {
  }
  
  template<typename TimeSteppingExplicitType>
  void TimeSteppingSchemeOdeReducedExplicit<TimeSteppingExplicitType>::
  initialize()
  {
    LOG(TRACE) << "TimeSteppingSchemeOdeReducedExplicit::initialize()";
    if (this->initialized_)
    {
      LOG(DEBUG) << "TimeSteppingSchemeOdeReducedExplicit is already initialized";
      return;
    }
    
    TimeSteppingSchemeOdeReduced<TimeSteppingExplicitType>::initialize(); 
    
    this->initialized_ = true;
  }
  
  template<typename TimeSteppingExplicitType>
  void TimeSteppingSchemeOdeReducedExplicit<TimeSteppingExplicitType>::
  evaluateTimesteppingRightHandSideExplicit(Vec &input, Vec &output, int timeStepNo, double currentTime)
  {
    this->fullTimestepping_.discretizableInTime().evaluateTimesteppingRightHandSideExplicit(input, output, timeStepNo, currentTime);
    
    // for implementing DEIM
    //ISCreateGeneral(comm,n, idx[],mode,&is);
    //VecGetSubVector(X,is,&Y) should it be called every time step?
  }
  
} //namespace
