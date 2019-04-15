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

namespace TimeSteppingScheme
{
/** This is the base class for all ode solvers.
 */
template<typename DiscretizableInTimeType>
class TimeSteppingSchemeOdeBaseDiscretizable:
public TimeSteppingSchemeOdeBase<typename DiscretizableInTimeType::FunctionSpace, DiscretizableInTimeType::nComponents()>
{
public:
  typedef DiscretizableInTimeType DiscretizableInTime_Type;
  typedef typename DiscretizableInTimeType::FunctionSpace FunctionSpace;
  
  //! constructor
  TimeSteppingSchemeOdeBaseDiscretizable(DihuContext context, std::string name);

  //! destructor
  virtual ~TimeSteppingSchemeOdeBaseDiscretizable() {}

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
template<int nStates, typename FunctionSpaceType>
class TimeSteppingSchemeOde<CellmlAdapter<nStates, FunctionSpaceType>>:
  public TimeSteppingSchemeOdeBaseDiscretizable<CellmlAdapter<nStates, FunctionSpaceType>>
{
public:
  //! use constructor of parent class
  using TimeSteppingSchemeOdeBaseDiscretizable<CellmlAdapter<nStates, FunctionSpaceType>>::TimeSteppingSchemeOdeBaseDiscretizable;
  
  //! initialize CellMLAdapter and get outputStateIndex and prefactor from CellMLAdapter to set them in data_
  virtual void initialize();
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


}  // namespace

#include "time_stepping_scheme/time_stepping_scheme_ode.tpp"
