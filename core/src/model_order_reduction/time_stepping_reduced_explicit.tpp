#include "model_order_reduction/time_stepping_reduced_explicit.h"

#include <Python.h>
#include<petscmat.h>
#include "utility/python_utility.h"
#include "utility/petsc_utility.h"


namespace ModelOrderReduction
{
  template<typename TimeSteppingExplicitType>
  TimeSteppingSchemeOdeReducedExplicit<TimeSteppingExplicitType>::
  TimeSteppingSchemeOdeReducedExplicit(DihuContext context):
  TimeSteppingSchemeOdeReduced<TimeSteppingExplicitType>(context,"ExplicitEulerReduced"), initialized_(false) 
  {
  }
  
  template<typename TimeSteppingExplicitType>
  void TimeSteppingSchemeOdeReducedExplicit<TimeSteppingExplicitType>::
  initialize()
  {
    if (initialized_)
      return;
    
    LOG(TRACE) << "TimeSteppingSchemeOdeReducedExplicit::initialize()";
    
    TimeSteppingSchemeOdeReduced<TimeSteppingExplicitType>::initialize(); 
    
    initialized_=true;
  }
  
  template<typename TimeSteppingExplicitType>
  void TimeSteppingSchemeOdeReducedExplicit<TimeSteppingExplicitType>::
  evaluateTimesteppingRightHandSideExplicit(Vec &input, Vec &output, int timeStepNo, double currentTime)
  {
    this->fullTimestepping_.discretizableInTime().evaluateTimesteppingRightHandSideExplicit(input, output, timeStepNo, currentTime);
    //ISCreateGeneral(comm,n, idx[],mode,&is);
    //VecGetSubVector(X,is,&Y) should it be called every time step?
  }
  
} //namespace
