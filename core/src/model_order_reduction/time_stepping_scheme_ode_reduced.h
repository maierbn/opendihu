#pragma once

#include "control/dihu_context.h"

#include "data_management/solution_vector_mapping.h"
#include "data_management/solution_vector_mapping.h"

#include "function_space/function_space.h"
#include "time_stepping_scheme/time_stepping_scheme.h"

#include "model_order_reduction/mor.h"

namespace ModelOrderReduction
{

template<typename TimeSteppingType>
class TimeSteppingSchemeOdeReduced :
  public MORBase<typename TimeSteppingType::FunctionSpace>,
  public TimeSteppingScheme::TimeSteppingScheme
{
public:
  typedef FieldVariable::FieldVariable<FunctionSpace::Generic,1> FieldVariableType;

  //! constructor
  TimeSteppingSchemeOdeReduced(DihuContext context);

  //! destructor
  virtual ~TimeSteppingSchemeOdeReduced(){};

  //! run simulation
  virtual void run();

  //! initialize timestepping member
  virtual void initialize();

  //!
  virtual void advanceTimeSpan()=0;

  //! return the Petsc solution vector
  std::shared_ptr<FieldVariableType> &solution();

  //! set the subset of ranks that will compute the work
  void setRankSubset(Partition::RankSubset rankSubset);

  //! reset state such that new initialization becomes necessary
  virtual void reset();

  ///! return whether the scheme has a specified mesh type and is not independent of the mesh type
  bool knowsMeshType();

protected:
  //! read initial values from settings and set field accordingly
  void setInitialValues();

  TimeSteppingType timestepping_;

private:
  bool initialized_;     ///< if initialize() was already called

};
  
} // namespace

#include "model_order_reduction/time_stepping_scheme_ode_reduced.tpp"
