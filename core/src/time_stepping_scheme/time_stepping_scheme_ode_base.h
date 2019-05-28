#pragma once

#include "control/dihu_context.h"
#include "data_management/time_stepping/time_stepping.h"
#include "interfaces/discretizable_in_time.h"
#include "time_stepping_scheme/time_stepping_scheme.h"
#include "data_management/data.h"
#include "data_management/time_stepping/time_stepping.h"
#include "cellml/03_cellml_adapter.h"

namespace TimeSteppingScheme
{


/**
  * Normal timestepping scheme to solve ODES, like Heun
  */
template<typename FunctionSpaceType, int nComponents>
class TimeSteppingSchemeOdeBase :
public TimeSteppingScheme
{
public:
  typedef FunctionSpaceType FunctionSpace;
  typedef Data::TimeStepping<FunctionSpaceType, nComponents> Data;   // type of Data object
  //typedef typename Data::TransferableSolutionDataType TransferableSolutionDataType;

  //! constructor
  TimeSteppingSchemeOdeBase(DihuContext context, std::string name);

  //! destructor
  virtual ~TimeSteppingSchemeOdeBase() {}

  //! run simulation
  virtual void run();

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the solution_vector_mapping class
  //virtual TransferableSolutionDataType getSolutionForTransfer();

  //! return the data object
  Data &data();

  //! output the given data for debugging
  //virtual std::string getString(TransferableSolutionDataType &data);

  //! initialize discretizableInTime
  virtual void initialize();

  //! set the subset of ranks that will compute the work
  void setRankSubset(Partition::RankSubset rankSubset);

  //! reset state such that new initialization becomes necessary
  virtual void reset();

protected:

  //! read initial values from settings and set field accordingly
  void setInitialValues();

  std::shared_ptr<Data> data_;     ///< data object that holds all PETSc vectors and matrices

  bool initialized_;     ///< if initialize() was already called

  double prefactor_;     ///< a factor with which the result is multiplied when the data is used in a splitting scheme

};

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

#include "time_stepping_scheme/time_stepping_scheme_ode_base.tpp"
#include "time_stepping_scheme/time_stepping_scheme_ode_transferable_solution_data.tpp"
