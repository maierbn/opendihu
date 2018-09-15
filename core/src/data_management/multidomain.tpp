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

template<typename FunctionSpaceType>
Multidomain<FunctionSpaceType>::
Multidomain(DihuContext context, int nCompartments) :
  Data<FunctionSpaceType>::Data(context), nCompartments_(nCompartments)
{

}


template<typename FunctionSpaceType>
void Multidomain<FunctionSpaceType>::
createPetscObjects()
{
  //this->boundaryConditionsRightHandSideSummand_ = this->functionSpace_->template createFieldVariable<nComponents>("boundaryConditionsRightHandSideSummand");
  this->fibreDirection_ = this->functionSpace_->template createFieldVariable<3>("fibreDirection");

}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> Multidomain<FunctionSpaceType>::
fibreDirection()
{
  return this->fibreDirection_;
}

template<typename FunctionSpaceType>
void Multidomain<FunctionSpaceType>::
print() // use override in stead of extending the parents' print output.This way "solution" is still in the end.
{
  if (!VLOG_IS_ON(4))
    return;

  VLOG(4) << *this->fibreDirection_;
}

} // namespace
