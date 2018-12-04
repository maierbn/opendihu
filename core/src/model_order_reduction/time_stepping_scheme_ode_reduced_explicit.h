#pragma once 

#include "control/dihu_context.h"
#include "model_order_reduction/time_stepping_scheme_ode_reduced.h"



namespace ModelOrderReduction
{
  template<typename TimeSteppingExplicitType>
  class TimeSteppingSchemeOdeReducedExplicit : 
  public TimeSteppingSchemeOdeReduced<TimeSteppingExplicitType>
  {
  public:
    
    //! constructor
    TimeSteppingSchemeOdeReducedExplicit(DihuContext context);
    
    //! destructor
    virtual ~TimeSteppingSchemeOdeReducedExplicit(){};
    
    //! run simulation
    void run();
    
    //! initialize timestepping member
    void initialize();
    
    //! advance the simulation by the time step
    void advanceTimeSpan();
    
    //! evaluates the right hand side function 
    virtual void evaluateTimesteppingRightHandSideExplicit(Vec &input, Vec &output, int timeStepNo, double currentTime);
    
  protected:
    
  private:
    bool initialized_;     ///< if initialize() was already called
  };
  
} // namespace

#include "model_order_reduction/time_stepping_scheme_ode_reduced_explicit.tpp"
