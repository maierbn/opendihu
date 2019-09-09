#pragma once

#include "control/dihu_context.h"
#include "data_management/time_stepping/time_stepping.h"
#include "interfaces/discretizable_in_time.h"
#include "time_stepping_scheme/time_stepping_scheme.h"
#include "data_management/data.h"
#include "data_management/time_stepping/time_stepping.h"
#include "cellml/03_cellml_adapter.h"
#include "spatial_discretization/boundary_conditions/dirichlet_boundary_conditions.h"
#include "time_stepping_scheme_ode_base.h"
//#include "time_stepping_scheme/time_stepping_scheme_ode_transferable_solution_data.h"

namespace TimeSteppingScheme
{
/** This is the base class for all ode solvers.
 */
template<typename DiscretizableInTimeType>
class TimeSteppingSchemeOdeBaseDiscretizable:
//  public TimeSteppingSchemeOdeOutputConnectorDataType<typename DiscretizableInTimeType::FunctionSpace, DiscretizableInTimeType::nComponents(), DiscretizableInTimeType>
  public TimeSteppingSchemeOdeBase<typename DiscretizableInTimeType::FunctionSpace, DiscretizableInTimeType::nComponents()>
{
public:
  typedef DiscretizableInTimeType DiscretizableInTime_Type;
  typedef typename DiscretizableInTimeType::FunctionSpace FunctionSpace;
  typedef typename DiscretizableInTimeType::OutputConnectorDataType OutputConnectorDataType;

  //using TimeSteppingSchemeOdeOutputConnectorDataType<typename DiscretizableInTimeType::FunctionSpace, DiscretizableInTimeType::nComponents(), DiscretizableInTimeType>::OutputConnectorDataType;

  //! constructor
  TimeSteppingSchemeOdeBaseDiscretizable(DihuContext context, std::string name);

  //! destructor
  virtual ~TimeSteppingSchemeOdeBaseDiscretizable() {}

  //! initialize discretizableInTime
  virtual void initialize();
  
  //! discretizable in time object
  DiscretizableInTimeType &discretizableInTime();
  
  //! Get the data that will be transferred in the operator splitting to the other term of the splitting.
  //! The transfer is done by the solution_vector_mapping class.
  //! The data is passed on from the DiscretizableInTime object.
  virtual OutputConnectorDataType getOutputConnectorData();

  //! set the subset of ranks that will compute the work
  void setRankSubset(Partition::RankSubset rankSubset);
  
  //! reset state such that new initialization becomes necessary
  virtual void reset();

  //! object that stores Dirichlet boundary condition values
  std::shared_ptr<
    SpatialDiscretization::DirichletBoundaryConditions<FunctionSpace,DiscretizableInTimeType::nComponents()>
  > dirichletBoundaryConditions();

protected:

  //! read initial values from settings and set field accordingly
  void setInitialValues();

  //int timeStepOutputInterval_;    ///< time step number and time is output every timeStepOutputInterval_ time steps
  DiscretizableInTimeType discretizableInTime_;    ///< the object to be discretized
  bool initialized_;     ///< if initialize() was already called

  std::shared_ptr<
    SpatialDiscretization::DirichletBoundaryConditions<FunctionSpace,DiscretizableInTimeType::nComponents()>
  > dirichletBoundaryConditions_;  ///< object that stores Dirichlet boundary condition values
};

template<typename DiscretizableInTimeType>
class TimeSteppingSchemeOde :
  public TimeSteppingSchemeOdeBaseDiscretizable<DiscretizableInTimeType>
  {
public:
  //! use constructor of parent class
  using TimeSteppingSchemeOdeBaseDiscretizable<DiscretizableInTimeType>::TimeSteppingSchemeOdeBaseDiscretizable;
  
  //using TimeSteppingSchemeOdeBaseDiscretizable<DiscretizableInTimeType>::initialize;
};

/**
 * Specialization for CellmlAdapter
 */
/*template<int nStates,int nIntermediates,typename FunctionSpaceType>
class TimeSteppingSchemeOde<CellmlAdapter<nStates, nIntermediates, FunctionSpaceType>> :
  public TimeSteppingSchemeOdeBaseDiscretizable<CellmlAdapter<nStates, nIntermediates, FunctionSpaceType>>
{
public:
  //! use constructor of parent class
  using TimeSteppingSchemeOdeBaseDiscretizable<CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>>::TimeSteppingSchemeOdeBaseDiscretizable;
  
  //! use data type for solution transfer
  typedef typename TimeSteppingSchemeOdeBaseDiscretizable<CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>>::OutputConnectorDataType OutputConnectorDataType;

  //! initialize CellMLAdapter and get outputStateIndex and prefactor from CellMLAdapter to set them in data_
  virtual void initialize();

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the solution_vector_mapping class
  virtual OutputConnectorDataType getOutputConnectorData();

  //! output the given data for debugging
  virtual std::string getString(OutputConnectorDataType &data);

protected:
  int outputIntermediateIndex_ = -1;   ///< component index of the intermediates field variable to use for solution transfer
};*/

}  // namespace

#include "time_stepping_scheme/time_stepping_scheme_ode.tpp"
#include "time_stepping_scheme/time_stepping_scheme_ode_cellml.tpp"
