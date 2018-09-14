#include "data_management/time_stepping_implicit.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <numeric>
#include <memory>

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "control/dihu_context.h"
#include "utility/petsc_utility.h"

namespace Data
{

template<typename FunctionSpaceType,int nComponents>
void TimeSteppingImplicit<FunctionSpaceType,nComponents>::
createPetscObjects()
{
  LOG(DEBUG) << "TimeSteppingImplicit::createPetscObjects(" <<nComponents << ")" << std::endl;
  TimeStepping<FunctionSpaceType,nComponents>::createPetscObjects();

  this->boundaryConditionsRightHandSideSummand_ = this->functionSpace_->template createFieldVariable<nComponents>("boundaryConditionsRightHandSideSummand");
  this->systemRightHandSide_ = this->functionSpace_->template createFieldVariable<nComponents>("systemRightHandSide");
}

template<typename FunctionSpaceType,int nComponents>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> TimeSteppingImplicit<FunctionSpaceType,nComponents>::
boundaryConditionsRightHandSideSummand()
{
  return this->boundaryConditionsRightHandSideSummand_;
}

template<typename FunctionSpaceType,int nComponents>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> TimeSteppingImplicit<FunctionSpaceType,nComponents>::
systemRightHandSide()
{
  return this->systemRightHandSide_;
}

template<typename FunctionSpaceType,int nComponents>
void TimeSteppingImplicit<FunctionSpaceType,nComponents>::
print() // use override in stead of extending the parents' print output.This way "solution" is still in the end.
{
  if (!VLOG_IS_ON(4))
    return;

  TimeStepping<FunctionSpaceType,nComponents>::print();

  VLOG(4) << *this->boundaryConditionsRightHandSideSummand_;
}

} // namespace
