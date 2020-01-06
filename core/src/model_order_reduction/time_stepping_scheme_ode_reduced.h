#pragma once

#include "control/dihu_context.h"

#include "function_space/function_space.h"
#include "time_stepping_scheme/02_time_stepping_scheme_ode.h"
#include "data_management/output_connector_data.h"
#include "data_management/time_stepping/time_stepping_reduced.h"
#include "model_order_reduction/model_order_reduction.h"

namespace ModelOrderReduction
{

  template<typename TimeSteppingType>
  class TimeSteppingSchemeOdeReduced :
    public MORBase<typename TimeSteppingType::FunctionSpace>,
    public ::TimeSteppingScheme::TimeSteppingSchemeOdeBase<::FunctionSpace::Generic,1>
  {
  public:
    typedef FieldVariable::FieldVariable<::FunctionSpace::Generic,1> FieldVariableType;  
    typedef ::FunctionSpace::Generic GenericFunctionSpace;
    typedef TimeSteppingType FullTimeSteppingType;
    typedef ::Data::OutputConnectorData<GenericFunctionSpace,1> OutputConnectorDataType;

    //! constructor
    TimeSteppingSchemeOdeReduced(DihuContext context,std::string name);

    //! destructor
    virtual ~TimeSteppingSchemeOdeReduced(){};
    
    //! run simulation
    virtual void run();

    //! initialize timestepping member
    virtual void initialize();

    //! set the subset of ranks that will compute the work
    void setRankSubset(Partition::RankSubset rankSubset);

    //! reset state such that new initialization becomes necessary
    //virtual void reset();

    //! full-order timestepping object
    TimeSteppingType fullTimestepping();

    //! get the data that will be transferred in the operator splitting to the other term of the splitting
    //! the transfer is done by the output_connector_data_transfer class
    std::shared_ptr<OutputConnectorDataType> getOutputConnectorData();

  protected:
    //! read initial values from settings and set field accordingly
    void setInitialValues();

    //! prepare the discretizableInTime object for the following call to getOutputConnectorData()
    virtual void prepareForGetOutputConnectorData() {}

    std::shared_ptr<GenericFunctionSpace> functionSpaceRed;
    std::shared_ptr<GenericFunctionSpace> functionSpaceRowsSnapshots;
    
    TimeSteppingType fullTimestepping_;
    std::shared_ptr<OutputConnectorDataType> outputConnectorData_;

    bool initialized_;     ///< if initialize() was already called

  };
  
} // namespace

#include "model_order_reduction/time_stepping_scheme_ode_reduced.tpp"
