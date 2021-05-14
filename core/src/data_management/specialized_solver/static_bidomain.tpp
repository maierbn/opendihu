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

  slotConnectorData_ = std::make_shared<SlotConnectorDataType>();
  slotConnectorData_->addFieldVariable(this->transmembranePotential_);

  // parse slot names for all slot connector data slots, only one slot here
  this->context_.getPythonConfig().getOptionVector("slotNames", slotConnectorData_->slotNames);

  // create additional field variables that appear as connector slots and can be connected to discretizableInTime_ and enclosing solvers
  int nAdditionalFieldVariables = this->context_.getPythonConfig().getOptionInt("nAdditionalFieldVariables", 0, PythonUtility::NonNegative);
  additionalFieldVariables_.resize(nAdditionalFieldVariables);

  for (int i = 0; i < nAdditionalFieldVariables; i++)
  {
    std::stringstream name;
    name << "additionalFieldVariable" << i;
    additionalFieldVariables_[i] = this->functionSpace_->template createFieldVariable<1>(name.str());

    slotConnectorData_->addFieldVariable2(additionalFieldVariables_[i]);
    LOG(DEBUG) << "  add field variable " << name.str();
  }
  
  // make sure that there are as many slot names as slots
  slotConnectorData_->slotNames.resize(slotConnectorData_->nSlots());
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
  this->jacobianConditionNumber_ = this->functionSpace_->template createFieldVariable<1>("jacobianConditionNumber");

  LOG(DEBUG) << "Vm field variable (" << this->transmembranePotential_ << ")";
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
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> StaticBidomain<FunctionSpaceType>::
jacobianConditionNumber()
{
  return this->jacobianConditionNumber_;
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
std::shared_ptr<typename StaticBidomain<FunctionSpaceType>::SlotConnectorDataType> StaticBidomain<FunctionSpaceType>::
getSlotConnectorData()
{
  return this->slotConnectorData_;
}

template<typename FunctionSpaceType>
typename StaticBidomain<FunctionSpaceType>::FieldVariablesForOutputWriter StaticBidomain<FunctionSpaceType>::
getFieldVariablesForOutputWriter()
{
  // these field variables will be written to output files
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> geometryField
    = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,3>>(this->functionSpace_->geometryField());

  // recover additional field variables from slotConnectorData_, they may have been changed by transfer
  assert(slotConnectorData_->variable2.size() >= additionalFieldVariables_.size());
  for (int i = 0; i < additionalFieldVariables_.size(); i++)
  {
    LOG(DEBUG) << " Data::StaticBidomain::getFieldVariablesForOutputWriter(), "
      << " get field variable " << slotConnectorData_->variable2[i].values << ", \"" << slotConnectorData_->variable2[i].values->name()
      << "\" for additionalFieldVariables_[" << i << "]";
    additionalFieldVariables_[i] = slotConnectorData_->variable2[i].values;
  }

  // these field variables will be written to output files
  return std::make_tuple(
    geometryField,
    this->fiberDirection_,
    extraCellularPotential_,
    transmembranePotential_,
    transmembraneFlow_,
    this->flowPotential_,
    jacobianConditionNumber_,
    additionalFieldVariables_
  );
}

} // namespace
