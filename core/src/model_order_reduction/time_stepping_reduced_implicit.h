#pragma once

#include "control/dihu_context.h"
#include "data_management/time_stepping/time_stepping_implicit.h"
#include "model_order_reduction/time_stepping_scheme_ode_reduced.h"

namespace ModelOrderReduction
{
  template<typename TimeSteppingImplicitType>
  class TimeSteppingSchemeOdeReducedImplicit : 
  public TimeSteppingSchemeOdeReduced<TimeSteppingImplicitType>,
  public ::TimeSteppingScheme::TimeSteppingImplicit<typename TimeSteppingImplicitType::DiscretizableInTime_Type>
  {
  public:
    typedef typename TimeSteppingSchemeOdeReduced<TimeSteppingImplicitType>::FunctionSpace FunctionSpace;
    
    //! constructor
    TimeSteppingSchemeOdeReducedImplicit(DihuContext context);
    
    //! destructor
    virtual ~TimeSteppingSchemeOdeReducedImplicit(){};
    
    //! advance simulation by the given time span [startTime_, endTime_] with given numberTimeSteps, data in solution is used, afterwards new data is in solution
    virtual void advanceTimeSpan() = 0;
    
    //! initialize timestepping member and set the system matrix
    void initialize();
    
  protected:
    
    //! precomputes the integration matrix for example A = (I-dtM^(-1)K) for the implicit euler scheme
    virtual void setSystemMatrix(double timeStepWidth) = 0;
    
    //! initialize the linear solve that is needed for the solution of the implicit timestepping system
    void initializeLinearSolver();
    
    //! solves the linear system of equations resulting from the Implicit Euler method time discretization
    void solveLinearSystem(Vec &input, Vec &output); 
    
    std::shared_ptr<Data::TimeSteppingImplicit<typename DiscretizableInTimeType::FunctionSpace, DiscretizableInTimeType::nComponents()>> dataImplicit_;  ///< a pointer to the data_ object but of type Data::TimeSteppingImplicit
    std::shared_ptr<Solver::Linear> linearSolver_;   ///< the linear solver used for solving the system
    std::shared_ptr<KSP> ksp_;     ///< the ksp object of the linear solver
    
  };
  
} // namespace

#include "model_order_reduction/time_stepping_reduced_implicit.tpp"
