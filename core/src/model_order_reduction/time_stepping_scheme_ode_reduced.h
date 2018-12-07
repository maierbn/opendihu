#pragma once

#include "control/dihu_context.h"

#include "function_space/function_space.h"
#include "time_stepping_scheme/time_stepping_scheme_ode.h"
#include "data_management/time_stepping/time_stepping_reduced.h"
#include "model_order_reduction/mor.h"

namespace ModelOrderReduction
{

  template<typename TimeSteppingType>
class TimeSteppingSchemeOdeReduced :
  public MORBase<typename TimeSteppingType::FunctionSpace>,
    public TimeSteppingScheme::TimeSteppingSchemeOdeBase<typename TimeSteppingType::DiscretizableInTime_Type>
{
public:
  //typedef ::FunctionSpace::Generic FunctionSpace;
  typedef FieldVariable::FieldVariable<::FunctionSpace::Generic,1> FieldVariableType;
  //typedef Data::TimeSteppingReduced<typename TimeSteppingType::FunctionSpace> DataReduced; //type of Data object  
  typedef FunctionSpace::Generic GenericFunctionSpace;
  
  //! constructor
  TimeSteppingSchemeOdeReduced(DihuContext context,std::string name);

  //! destructor
  virtual ~TimeSteppingSchemeOdeReduced(){};
  
  //! run simulation
  virtual void run();

  //! initialize timestepping member
  virtual void initialize();

  //!
  virtual void advanceTimeSpan()=0;
  
  //! data object for model order reduction
  //DataReduced &redData();

  //! return the Petsc solution vector
  std::shared_ptr<FieldVariableType> &solution();
  
  //! return the data object
  //DataReduced &redData();

  //! set the subset of ranks that will compute the work
  void setRankSubset(Partition::RankSubset rankSubset);

  //! reset state such that new initialization becomes necessary
  //virtual void reset();

  //! return whether the scheme has a specified mesh type and is not independent of the mesh type
  bool knowsMeshType();

protected:
  //! read initial values from settings and set field accordingly
  void setInitialValues();
  
  //std::shared_ptr<DataReduced> redData_;     ///< data object that holds all PETSc vectors and matrices
  
  //std::shared_ptr<DataReduced> redData_;
  
  std::shared_ptr<GenericFunctionSpace> functionSpaceRed;
  std::shared_ptr<GenericFunctionSpace> functionSpaceRowsSnapshots;
  
  TimeSteppingType fullTimestepping_;

private:
  bool initialized_;     ///< if initialize() was already called

};
  
} // namespace

#include "model_order_reduction/time_stepping_scheme_ode_reduced.tpp"
