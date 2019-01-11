#pragma once 

#include "control/dihu_context.h"
#include "function_space/function_space.h"
#include "model_order_reduction/time_stepping_scheme_ode_reduced.h"

namespace ModelOrderReduction
{
  template<typename TimeSteppingExplicitType>
  class TimeSteppingSchemeOdeReducedExplicit : 
  public TimeSteppingSchemeOdeReduced<TimeSteppingExplicitType>
  {
  public:
    typedef typename TimeSteppingSchemeOdeReduced<TimeSteppingExplicitType>::FunctionSpace FunctionSpace;
    
    //! constructor
    TimeSteppingSchemeOdeReducedExplicit(DihuContext context);
    
    //! destructor
    virtual ~TimeSteppingSchemeOdeReducedExplicit(){};
    
    //! initialize timestepping member
    void initialize();
    
    //! evaluates the right hand side function 
    void evaluateTimesteppingRightHandSideExplicit(Vec &input, Vec &output, int timeStepNo, double currentTime);
    
  protected:    
    
  };
  
} // namespace

#include "model_order_reduction/time_stepping_reduced_explicit.tpp"
