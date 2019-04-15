#pragma once

#include "control/dihu_context.h"
#include "data_management/time_stepping/time_stepping.h"
#include "interfaces/discretizable_in_time.h"
#include "time_stepping_scheme/time_stepping_scheme.h"
#include "time_stepping_scheme/time_stepping_scheme_ode.h"
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

protected:

  //! set the dofs in the solution vector to the given boundary conditions
  void applyBoundaryConditions();
};
}  // namespace

#include "time_stepping_scheme/time_stepping_explicit.tpp"
