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
class TimeSteppingExplicit:
  public TimeSteppingSchemeOde<DiscretizableInTimeType>
{
public:
  //! constructor
  using TimeSteppingSchemeOde<DiscretizableInTimeType>::TimeSteppingSchemeOde;

  //! initialize discretizableInTime
  virtual void initialize();

protected:

  //! set the dofs in the solution vector to the given boundary conditions
  void applyBoundaryConditions();

  std::shared_ptr<SpatialDiscretization::DirichletBoundaryConditions<typename DiscretizableInTimeType::FunctionSpace, DiscretizableInTimeType::nComponents()>> dirichletBoundaryConditions_; ///< the dirichlet boundary conditions object that parses, stores and applies the dirichlet boundary conditions
};
}  // namespace

#include "time_stepping_scheme/time_stepping_explicit.tpp"
