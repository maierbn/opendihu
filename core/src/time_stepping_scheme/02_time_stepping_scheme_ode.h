#pragma once

#include <Python.h>  // has to be the first included header

#include "control/dihu_context.h"
#include "data_management/time_stepping/time_stepping.h"
#include "interfaces/discretizable_in_time.h"
#include "time_stepping_scheme/00_time_stepping_scheme.h"
#include "data_management/data.h"
#include "data_management/time_stepping/time_stepping.h"
#include "cellml/03_cellml_adapter.h"
#include "spatial_discretization/dirichlet_boundary_conditions/01_dirichlet_boundary_conditions.h"
#include "time_stepping_scheme/01_time_stepping_scheme_ode_base.h"
//#include "time_stepping_scheme/time_stepping_scheme_ode_transferable_solution_data.h"

namespace TimeSteppingScheme
{
/** This is the base class for all ode solvers.
 */
template<typename DiscretizableInTimeType>
class TimeSteppingSchemeOdeBaseDiscretizable:
  public TimeSteppingSchemeOdeBase<typename DiscretizableInTimeType::FunctionSpace, DiscretizableInTimeType::nComponents()>
{
public:
  typedef DiscretizableInTimeType DiscretizableInTime;
  typedef typename DiscretizableInTimeType::FunctionSpace FunctionSpace;
  //typedef typename DiscretizableInTimeType::SlotConnectorDataType SlotConnectorDataType;

  //using TimeSteppingSchemeOdeSlotConnectorDataType<typename DiscretizableInTimeType::FunctionSpace, DiscretizableInTimeType::nComponents(), DiscretizableInTimeType>::SlotConnectorDataType;

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

  //! object that stores Dirichlet boundary condition values
  std::shared_ptr<
    SpatialDiscretization::DirichletBoundaryConditions<FunctionSpace,DiscretizableInTimeType::nComponents()>
  > dirichletBoundaryConditions();

protected:

  //! read initial values from settings and set field accordingly
  void setInitialValues();

  //! prepare the discretizableInTime object for the following call to getSlotConnectorData()
  virtual void prepareForGetSlotConnectorData() override;

  //int timeStepOutputInterval_;    //< time step number and time is output every timeStepOutputInterval_ time steps
  DiscretizableInTimeType discretizableInTime_;    //< the object to be discretized
  bool initialized_;     //< if initialize() was already called

  std::shared_ptr<
    SpatialDiscretization::DirichletBoundaryConditions<FunctionSpace,DiscretizableInTimeType::nComponents()>
  > dirichletBoundaryConditions_;  //< object that stores Dirichlet boundary condition values
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
}  // namespace

#include "time_stepping_scheme/02_time_stepping_scheme_ode.tpp"
