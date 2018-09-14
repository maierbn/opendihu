#include "time_stepping_scheme/time_stepping_explicit.h"

#include <Python.h>  // has to be the first included header
#include <vector>

#include "utility/python_utility.h"
#include "spatial_discretization/dirichlet_boundary_conditions.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTimeType>
void TimeSteppingExplicit<DiscretizableInTimeType>::
initialize()
{
  TimeSteppingSchemeOde<DiscretizableInTimeType>::initialize();

  // parse boundary conditions, needs functionSpace set
  dirichletBoundaryConditions_ = std::make_shared<::SpatialDiscretization::DirichletBoundaryConditions<typename DiscretizableInTimeType::FunctionSpace, DiscretizableInTimeType::nComponents()>>();
  dirichletBoundaryConditions_->initialize(this->specificSettings_, this->data_->functionSpace());
}

template<typename DiscretizableInTimeType>
void TimeSteppingExplicit<DiscretizableInTimeType>::
applyBoundaryConditions()
{
  dirichletBoundaryConditions_->applyInVector(this->data_->solution());
}

} // namespace
