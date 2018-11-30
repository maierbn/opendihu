#pragma once

#include "control/dihu_context.h"
#include "model_order_reduction/time_stepping_scheme_ode_reduced.h"

namespace ModelOrderReduction
{
  template<typename TimeSteppingImplicitType>
  class TimeSteppingSchemeOdeReducedImplicit : 
  public TimeSteppingSchemeOdeReduced<TimeSteppingImplicitType>
  {
  public:
    typedef typename TimeSteppingSchemeOdeReduced<TimeSteppingImplicitType>::FunctionSpace FunctionSpace;
    
    //! constructor
    TimeSteppingSchemeOdeReducedImplicit(DihuContext context);
    
    //! destructor
    virtual ~TimeSteppingSchemeOdeReducedImplicit(){};
    
    //! run simulation
    void run();
    
    //! initialize timestepping member
    void initialize();
    
    //! 
    void advanceTimeSpan();
    
  protected:
    
  private:
    bool initialized_;     ///< if initialize() was already called
  };
  
} // namespace

#include "model_order_reduction/time_stepping_reduced_implicit.tpp"
