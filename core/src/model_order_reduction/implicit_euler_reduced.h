#pragma once

#include "control/dihu_context.h"
#include "model_order_reduction/time_stepping_scheme_reduced_implicit.h"

namespace ModelOrderReduction
{
  template<typename TimeSteppingImplicitType>
  class  : ImplicitEulerReduced
  public TimeSteppingSchemeOdeReducedImplicit<TimeSteppingImplicitType>
  {
  public:
    typedef typename TimeSteppingSchemeOdeReduced<TimeSteppingImplicitType>::FunctionSpace FunctionSpace;
    
    //! constructor
    ImplicitEulerReduced(DihuContext context);
    
    //! destructor
    virtual ~ImplicitEulerReduced(){};
    
    //! run simulation
    void run();
    
    //! 
    void advanceTimeSpan();
    
  protected:
    
  };
  
} // namespace

#include "model_order_reduction/implicit_euler_reduced.tpp"
