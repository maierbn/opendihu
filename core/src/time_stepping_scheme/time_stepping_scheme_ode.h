#pragma once

#include "control/dihu_context.h"
#include "data_management/solution_vector_mapping.h"
#include "data_management/time_stepping.h"
#include "discretizable_in_time/discretizable_in_time.h"
#include "time_stepping_scheme/time_stepping_scheme.h"
#include "data_management/data.h"

namespace TimeSteppingScheme
{

/** This is the base class for all ode solvers.
 */
template<typename DiscretizableInTimeType>
class TimeSteppingSchemeOde : public TimeSteppingScheme
{
public:
  typedef typename DiscretizableInTimeType::FunctionSpace FunctionSpace;

  typedef Data::TimeStepping<typename DiscretizableInTimeType::FunctionSpace, DiscretizableInTimeType::nComponents()> Data;   // type of Data object

  //! constructor
  TimeSteppingSchemeOde(DihuContext context, const std::string name);

  //! destructor
  virtual ~TimeSteppingSchemeOde() {}

  //! run simulation
  virtual void run();

  //! get the solutionVectorMapping object that stores information about which values of the solution should be used for further computation and how they can be retrieved
  SolutionVectorMapping &solutionVectorMapping();

  //! return the Petsc solution vector
  typename Data::FieldVariableType &solution();

  //! return the data object
  Data &data();

  //! initialize discretizableInTime
  virtual void initialize();
  
  //! set the subset of ranks that will compute the work
  void setRankSubset(Partition::RankSubset rankSubset);
  
  //! reset state such that new initialization becomes necessary
  virtual void reset();

  //! return whether the underlying discretizableInTime object has a specified mesh type and is not independent of the mesh type
  bool knowsMeshType();

protected:

  //! read initial values from settings and set field accordingly
  void setInitialValues();

  std::shared_ptr<Data> data_;     ///< data object that holds all PETSc vectors and matrices

  int timeStepOutputInterval_;    ///< time step number and time is output every timeStepOutputInterval_ time steps
  DiscretizableInTimeType discretizableInTime_;    ///< the object to be discretized
  bool initialized_;     ///< if initialize() was already called
};
}  // namespace

#include "time_stepping_scheme/time_stepping_scheme_ode.tpp"
