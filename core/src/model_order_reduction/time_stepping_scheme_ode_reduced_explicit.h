#include "control/dihu_context.h"
#include "data_management/solution_vector_mapping.h"
#include "model_order_reduction/mor.h"
#include "time_stepping_scheme/time_stepping_scheme.h"



namespace ModelOrderReduction
{
  template<typename TimeSteppingExplicitType>
  class TimeSteppingSchemeOdeReducedExplicit <TimeSteppingExplicitType> : 
    public MORBase, 
    public TimeSteppingScheme::TimeSteppingScheme
  {
  public:
    typdef typename FieldVariable::FieldVariable<FunctionSpace::Generic,1> FieldVariableType;
    
    //! constructor
    TimeSteppingSchemeOdeReducedExplicit(DihuContext context);
    
    //! destructor
    virtual ~TimeSteppingSchemeOdeReducedExplicit(){};
    
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

#include "model_order_reduction/time_stepping_scheme_ode_reduced_explicit.tpp"