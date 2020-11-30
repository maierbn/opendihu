#pragma once

#include "control/dihu_context.h"
#include "model_order_reduction/time_stepping_reduced_implicit.h"

namespace ModelOrderReduction
{
  template<typename TimeSteppingImplicitType>
  class  ImplicitEulerReduced:
  public TimeSteppingSchemeOdeReducedImplicit<TimeSteppingImplicitType>
  {
  public:
    
    //! constructor
    ImplicitEulerReduced(DihuContext context);
    
    //! destructor
    virtual ~ImplicitEulerReduced(){};
    
    //! run simulation
    void run();
    
    //! 
    void advanceTimeSpan(bool withOutputWritersEnabled = true);
    
  protected:
    
  };
  
} // namespace

#include "model_order_reduction/implicit_euler_reduced.tpp"
