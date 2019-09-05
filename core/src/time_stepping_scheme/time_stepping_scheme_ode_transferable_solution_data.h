#pragma once

#include "time_stepping_scheme/time_stepping_scheme_ode_base.h"

namespace TimeSteppingScheme
{

template<typename FunctionSpaceType, int nComponents, typename DiscretizableInTimeType=FunctionSpaceType>
class TimeSteppingSchemeOdeOutputConnectorDataType :
  public TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>
{
public:
  using TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::TimeSteppingSchemeOdeBase;

  typedef typename Data::TimeStepping<FunctionSpaceType, nComponents> Data;   // type of Data object
  typedef typename Data::OutputConnectorDataType OutputConnectorDataType;

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the solution_vector_mapping class
  virtual OutputConnectorDataType getOutputConnectorData();

  //! output the given data for debugging
  virtual std::string getString(OutputConnectorDataType &data);

};

/** Partial specialization for CellML adapter, which transfers different data (states and intermediates)
 */
template<typename FunctionSpaceType,int nComponents,int nStates,int nIntermediates>
class TimeSteppingSchemeOdeOutputConnectorDataType<FunctionSpaceType, nComponents, CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>> :
public TimeSteppingSchemeOdeBase<FunctionSpaceType,nComponents>
{
public:
  using TimeSteppingSchemeOdeBase<FunctionSpaceType,nComponents>::TimeSteppingSchemeOdeBase;

  typedef typename Data::TimeStepping<FunctionSpaceType,nComponents> DataType;   // type of Data object
  typedef std::pair<
    typename DataType::OutputConnectorDataType,    // std::tuple<std::shared_ptr<FieldVariableType>,int,double> OutputConnectorDataType;  // <field variable, output component no., prefactor>
    std::tuple<std::shared_ptr<typename CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::FieldVariableIntermediates>,int>  // intermediates, output component no., -1 if unused
  > OutputConnectorDataType;

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the solution_vector_mapping class
  virtual OutputConnectorDataType getOutputConnectorData() = 0;

  //! output the given data for debugging
  virtual std::string getString(OutputConnectorDataType &data) = 0;
};

}  // namespace

#include "time_stepping_scheme/time_stepping_scheme_ode_transferable_solution_data.tpp"
