#include "data_management/specialized_solver/static_bidomain.h"

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
StaticBidomain<FunctionSpaceType>::
StaticBidomain(DihuContext context) :
  Data<FunctionSpaceType>::Data(context)
{
}

template<typename FunctionSpaceType>
void StaticBidomain<FunctionSpaceType>::
initialize()
{
  // call initialize of base class
  Data<FunctionSpaceType>::initialize();
}

template<typename FunctionSpaceType>
void StaticBidomain<FunctionSpaceType>::
createPetscObjects()
{
  LOG(DEBUG) << "StaticBidomain::createPetscObject";

  assert(this->functionSpace_);

  this->transmembraneFlow_ = this->functionSpace_->template createFieldVariable<1>("transmembraneFlow");
  this->transmembranePotential_ = this->functionSpace_->template createFieldVariable<1>("Vm");
  this->flowPotential_ = this->functionSpace_->template createFieldVariable<1>("flowPotential");
  this->fiberDirection_ = this->functionSpace_->template createFieldVariable<3>("fiberDirection");
  this->extraCellularPotential_ = this->functionSpace_->template createFieldVariable<1>("phi_e");
  this->zero_ = this->functionSpace_->template createFieldVariable<1>("zero");
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> StaticBidomain<FunctionSpaceType>::
fiberDirection()
{
  return this->fiberDirection_;
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> StaticBidomain<FunctionSpaceType>::
flowPotential()
{
  return this->flowPotential_;
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> StaticBidomain<FunctionSpaceType>::
extraCellularPotential()
{
  return this->extraCellularPotential_;
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> StaticBidomain<FunctionSpaceType>::
transmembranePotential()
{
  return this->transmembranePotential_;
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> StaticBidomain<FunctionSpaceType>::
transmembraneFlow()
{
  return this->transmembraneFlow_;
}

template<typename FunctionSpaceType>
Mat &StaticBidomain<FunctionSpaceType>::
rhsMatrix()
{
  return this->rhsMatrix_;
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> StaticBidomain<FunctionSpaceType>::
zero()
{
  return this->zero_;
}

template<typename FunctionSpaceType>
void StaticBidomain<FunctionSpaceType>::
print() // use override in stead of extending the parents' print output.This way "solution" is still in the end.
{
  if (!VLOG_IS_ON(4))
    return;

  VLOG(4) << *this->fiberDirection_;
}

template<typename FunctionSpaceType>
typename StaticBidomain<FunctionSpaceType>::FieldVariablesForOutputWriter StaticBidomain<FunctionSpaceType>::
getFieldVariablesForOutputWriter()
{
  // these field variables will be written to output files
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> geometryField
    = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,3>>(this->functionSpace_->geometryField());

  return std::make_tuple(
    geometryField,
    this->fiberDirection_,
    this->flowPotential_,
    extraCellularPotential_,
    transmembranePotential_,
    transmembraneFlow_
  );
}

} // namespace
