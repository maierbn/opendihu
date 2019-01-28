#pragma once

#include "control/dihu_context.h"
#include "data_management/time_stepping/time_stepping.h"
#include "interfaces/discretizable_in_time.h"
#include "time_stepping_scheme/time_stepping_scheme.h"
#include "data_management/data.h"
#include "cellml/03_cellml_adapter.h"

namespace TimeSteppingScheme
{

/** This is the base class for all ode solvers.
 */
template<typename DiscretizableInTimeType>
class TimeSteppingSchemeOdeBase :
  public TimeSteppingScheme
{
public:
  typedef typename DiscretizableInTimeType::FunctionSpace FunctionSpace;
  typedef Data::TimeStepping<typename DiscretizableInTimeType::FunctionSpace, DiscretizableInTimeType::nComponents()> Data;   // type of Data object
  typedef typename Data::TransferableSolutionDataType TransferableSolutionDataType;

  //! constructor
  TimeSteppingSchemeOdeBase(DihuContext context, std::string name);

  //! destructor
  virtual ~TimeSteppingSchemeOdeBase() {}

  //! run simulation
  virtual void run();

  //! get the data that will be transferred in the operator splitting to the other term of the splitting
  //! the transfer is done by the solution_vector_mapping class
  TransferableSolutionDataType getSolutionForTransfer();

  //! return the data object
  Data &data();

  //! initialize discretizableInTime
  virtual void initialize();
  
  //! discretizable in time object
  DiscretizableInTimeType &discretizableInTime();
  
  //! set the subset of ranks that will compute the work
  void setRankSubset(Partition::RankSubset rankSubset);
  
  //! reset state such that new initialization becomes necessary
  virtual void reset();

  //! return whether the underlying discretizableInTime object has a specified mesh type and is not independent of the mesh type
  bool knowsMeshType();

  //! output the given data for debugging
  virtual std::string getString(TransferableSolutionDataType &data);

protected:

  //! read initial values from settings and set field accordingly
  void setInitialValues();

  std::shared_ptr<Data> data_;     ///< data object that holds all PETSc vectors and matrices

  //int timeStepOutputInterval_;    ///< time step number and time is output every timeStepOutputInterval_ time steps
  DiscretizableInTimeType discretizableInTime_;    ///< the object to be discretized
  bool initialized_;     ///< if initialize() was already called

  double prefactor_;     ///< a factor with which the result is multiplied when the data is used in a splitting scheme

  std::shared_ptr<
    SpatialDiscretization::DirichletBoundaryConditions<FunctionSpace,DiscretizableInTimeType::nComponents()>
  > dirichletBoundaryConditions_;  ///< object that stores Dirichlet boundary condition values
};

template<typename DiscretizableInTimeType>
class TimeSteppingSchemeOde :
  public TimeSteppingSchemeOdeBase<DiscretizableInTimeType>
{
public:
  //! use constructor of parent class
  using TimeSteppingSchemeOdeBase<DiscretizableInTimeType>::TimeSteppingSchemeOdeBase;

};

/**
 * Specialization for CellmlAdapter
 */
template<int nStates, typename FunctionSpaceType>
class TimeSteppingSchemeOde<CellmlAdapter<nStates, FunctionSpaceType>> :
  public TimeSteppingSchemeOdeBase<CellmlAdapter<nStates, FunctionSpaceType>>
{
public:
  //! use constructor of parent class
  using TimeSteppingSchemeOdeBase<CellmlAdapter<nStates, FunctionSpaceType>>::TimeSteppingSchemeOdeBase;

  //! initialize CellMLAdapter and get outputStateIndex and prefactor from CellMLAdapter to set them in data_
  virtual void initialize();
};


}  // namespace

#include "time_stepping_scheme/time_stepping_scheme_ode.tpp"
