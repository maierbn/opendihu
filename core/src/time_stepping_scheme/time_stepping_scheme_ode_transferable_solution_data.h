#pragma once

#include "time_stepping_scheme/time_stepping_scheme_ode_base.h"

namespace TimeSteppingScheme
{

template<typename FunctionSpaceType, int nComponents, typename DiscretizableInTimeType=FunctionSpaceType>
class TimeSteppingSchemeOdeTransferableSolutionData :
  public TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>
{
public:
  using TimeSteppingSchemeOdeBase<FunctionSpaceType, nComponents>::TimeSteppingSchemeOdeBase;

  typedef Data::TimeStepping<FunctionSpaceType, nComponents> Data;   // type of Data object
  typedef typename Data::TransferableSolutionDataType TransferableSolutionDataType;

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the solution_vector_mapping class
  virtual TransferableSolutionDataType getSolutionForTransfer();

  //! output the given data for debugging
  virtual std::string getString(TransferableSolutionDataType &data);

};

/** Partial specialization for CellML adapter, which transfers different data (states and intermediates)
 */
template<typename FunctionSpaceType,int nComponents,int nStates,int nIntermediates>
class TimeSteppingSchemeOdeTransferableSolutionData<FunctionSpaceType, nComponents, CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>> :
public TimeSteppingSchemeOdeBase<FunctionSpaceType,nComponents>
{
public:
  using TimeSteppingSchemeOdeBase<FunctionSpaceType,nComponents>::TimeSteppingSchemeOdeBase;

  typedef Data::TimeStepping<FunctionSpaceType,nComponents> DataType;   // type of Data object
  typedef std::pair<
    typename DataType::TransferableSolutionDataType,    // std::tuple<std::shared_ptr<FieldVariableType>,int,double> TransferableSolutionDataType;  // <field variable, output component no., prefactor>
    std::tuple<std::shared_ptr<typename CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::FieldVariableTypeIntermediates>,int>  // intermediates, output component no., -1 if unused
  > TransferableSolutionDataType;

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the solution_vector_mapping class
  virtual TransferableSolutionDataType getSolutionForTransfer() = 0;

  //! output the given data for debugging
  virtual std::string getString(TransferableSolutionDataType &data) = 0;
};

}  // namespace

#include "time_stepping_scheme/time_stepping_scheme_ode_transferable_solution_data.tpp"
