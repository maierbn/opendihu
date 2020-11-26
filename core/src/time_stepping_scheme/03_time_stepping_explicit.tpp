#include "time_stepping_scheme/03_time_stepping_explicit.h"

#include <Python.h>  // has to be the first included header
#include <vector>

#include "utility/python_utility.h"

namespace TimeSteppingScheme
{

template<typename DiscretizableInTimeType>
void TimeSteppingExplicit<DiscretizableInTimeType>::
applyBoundaryConditions()
{
  this->dirichletBoundaryConditions_->applyInVector(this->data_->solution());
}

} // namespace
