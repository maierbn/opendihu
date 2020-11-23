#pragma once 

#include "control/dihu_context.h"
#include "function_space/function_space.h"
#include "model_order_reduction/time_stepping_reduced_explicit.h"

namespace ModelOrderReduction
{
  template<typename TimeSteppingExplicitType>
  class  ExplicitEulerReduced: 
  public TimeSteppingSchemeOdeReducedExplicit<TimeSteppingExplicitType>
  {
  public:
    
    //! constructor
    ExplicitEulerReduced(DihuContext context);
    
    //! destructor
    virtual ~ExplicitEulerReduced(){};
    
    //! run simulation
    void run();
    
    //! advance the simulation by the time step
    void advanceTimeSpan(bool withOutputWritersEnabled = true);
    
  protected:
  
  };
  
} // namespace

#include "model_order_reduction/explicit_euler_reduced.tpp"
