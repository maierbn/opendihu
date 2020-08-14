#include "data_management/specialized_solver/multidomain.h"

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
Multidomain(DihuContext context) :
  Data<FunctionSpaceType>::Data(context)
{
}

template<typename FunctionSpaceType>
void Multidomain<FunctionSpaceType>::
initialize(int nCompartments)
{
  if (this->initialized_)
    return;
    
  nCompartments_ = nCompartments;

  // call initialize of base class
  Data<FunctionSpaceType>::initialize();

  // initialize slot connector data
  std::get<0>(*slotConnectorData_) = std::make_shared<std::vector<std::shared_ptr<SlotConnectorData<FunctionSpaceType,1>>>>();
  std::get<0>(*slotConnectorData_)->resize(nCompartments_);

  std::string activeStressTotalSlotName;

  // set the two slots V_mk^(i) and V_mk^(i+1) for all compartments
  for (int compartmentNo = 0; compartmentNo < nCompartments_; compartmentNo++)
  {
    std::get<0>(*slotConnectorData_)->at(compartmentNo) = std::make_shared<SlotConnectorData<FunctionSpaceType,1>>();
    std::get<0>(*slotConnectorData_)->at(compartmentNo)->addFieldVariable(transmembranePotential_[compartmentNo]);          // V_mk^(i)
    std::get<0>(*slotConnectorData_)->at(compartmentNo)->addFieldVariable(transmembranePotentialSolution_[compartmentNo]);  // V_mk^(i+1)

    std::get<0>(*slotConnectorData_)->at(compartmentNo)->addFieldVariable(activeStress_[compartmentNo]);  // activeStress_k

    // parse slot names of the field variables
    this->context_.getPythonConfig().getOptionVector("slotNames", std::get<0>(*slotConnectorData_)->at(compartmentNo)->slotNames);

    int nSlots = std::get<0>(*slotConnectorData_)->at(compartmentNo)->nSlots();
    if (std::get<0>(*slotConnectorData_)->at(compartmentNo)->slotNames.size() >= nSlots+1)
    {
      activeStressTotalSlotName = std::get<0>(*slotConnectorData_)->at(compartmentNo)->slotNames[nSlots];
    }

    // make sure that there are as many slot names as slots
    std::get<0>(*slotConnectorData_)->at(compartmentNo)->slotNames.resize(nSlots);
  }

  // active stress total
  std::get<1>(*slotConnectorData_) = std::make_shared<SlotConnectorData<FunctionSpaceType,1>>();
  std::get<1>(*slotConnectorData_)->addFieldVariable(activeStressTotal_);  // activeStressTotal

  // parse slot names of the field variables
  this->context_.getPythonConfig().getOptionVector("slotNames", std::get<1>(*slotConnectorData_)->slotNames);

  LOG(DEBUG) << "add slot for activeStressTotal with slotName \"" << activeStressTotalSlotName << "\"";
  std::get<1>(*slotConnectorData_)->slotNames.push_back(activeStressTotalSlotName);
}

template<typename FunctionSpaceType>
void Multidomain<FunctionSpaceType>::
createPetscObjects()
{
  LOG(DEBUG) << "Multidomain::createPetscObject for " << nCompartments_ << " compartments.";

  // create field variables that have one for every compartment
  transmembranePotentialSolution_.reserve(nCompartments_);
  transmembranePotential_.reserve(nCompartments_);
  compartmentRelativeFactor_.reserve(nCompartments_);
  activeStress_.reserve(nCompartments_);

  assert(this->functionSpace_);

  for (int k = 0; k < nCompartments_; k++)
  {
    std::stringstream transmembranePotentialSolutionName;
    transmembranePotentialSolutionName << "Vm^(i+1)_" << k;
    this->transmembranePotentialSolution_.push_back(this->functionSpace_->template createFieldVariable<1>(transmembranePotentialSolutionName.str()));

    std::stringstream transmembranePotentialName;
    transmembranePotentialName << "Vm^(i)_" << k;
    this->transmembranePotential_.push_back(this->functionSpace_->template createFieldVariable<1>(transmembranePotentialName.str()));

    std::stringstream compartmentRelativeFactorName;
    compartmentRelativeFactorName << "f_r_" << k;
    this->compartmentRelativeFactor_.push_back(this->functionSpace_->template createFieldVariable<1>(compartmentRelativeFactorName.str()));

    std::stringstream activeStressName;
    activeStressName << "active_stress_" << k;
    this->activeStress_.push_back(this->functionSpace_->template createFieldVariable<1>(activeStressName.str()));
  }

  this->flowPotential_ = this->functionSpace_->template createFieldVariable<1>("flowPotential");
  this->fiberDirection_ = this->functionSpace_->template createFieldVariable<3>("fiberDirection");
  this->extraCellularPotential_ = this->functionSpace_->template createFieldVariable<1>("phi_e");
  this->zero_ = this->functionSpace_->template createFieldVariable<1>("zero");
  this->relativeFactorTotal_ = this->functionSpace_->template createFieldVariable<1>("Σf_r");
  this->activeStressTotal_ = this->functionSpace_->template createFieldVariable<1>("activeStressTotal");
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> Multidomain<FunctionSpaceType>::
fiberDirection()
{
  return this->fiberDirection_;
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> Multidomain<FunctionSpaceType>::
flowPotential()
{
  return this->flowPotential_;
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> Multidomain<FunctionSpaceType>::
extraCellularPotential()
{
  return this->extraCellularPotential_;
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> Multidomain<FunctionSpaceType>::
transmembranePotential(int compartmentNo)
{
  assert(compartmentNo >= 0 && compartmentNo < nCompartments_);
  return this->transmembranePotential_[compartmentNo];
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> Multidomain<FunctionSpaceType>::
transmembranePotentialSolution(int compartmentNo)
{
  assert(compartmentNo >= 0 && compartmentNo < nCompartments_);
  return this->transmembranePotentialSolution_[compartmentNo];
}

template<typename FunctionSpaceType>
std::vector<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>>> Multidomain<FunctionSpaceType>::
transmembranePotential()
{
  return this->transmembranePotential_;
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> Multidomain<FunctionSpaceType>::
compartmentRelativeFactor(int compartmentNo)
{
  assert(compartmentNo >= 0 && compartmentNo < nCompartments_);
  return this->compartmentRelativeFactor_[compartmentNo];
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> Multidomain<FunctionSpaceType>::
activeStress(int compartmentNo)
{
  assert(compartmentNo >= 0 && compartmentNo < nCompartments_);
  return this->activeStress_[compartmentNo];
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> Multidomain<FunctionSpaceType>::
relativeFactorTotal()
{
  return this->relativeFactorTotal_;
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> Multidomain<FunctionSpaceType>::
activeStressTotal()
{
  return this->activeStressTotal_;
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> Multidomain<FunctionSpaceType>::
zero()
{
  return this->zero_;
}

template<typename FunctionSpaceType>
std::shared_ptr<typename Multidomain<FunctionSpaceType>::SlotConnectorDataType> Multidomain<FunctionSpaceType>::
getSlotConnectorData()
{
  return slotConnectorData_;
}

template<typename FunctionSpaceType>
void Multidomain<FunctionSpaceType>::
print() // use override in stead of extending the parents' print output.This way "solution" is still in the end.
{
  if (!VLOG_IS_ON(4))
    return;

  VLOG(4) << *this->fiberDirection_;
}

template<typename FunctionSpaceType>
typename Multidomain<FunctionSpaceType>::FieldVariablesForOutputWriter Multidomain<FunctionSpaceType>::
getFieldVariablesForOutputWriter()
{
  // these field variables will be written to output files
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> geometryField
    = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,3>>(this->functionSpace_->geometryField());

  std::vector<std::shared_ptr<FieldVariableType>> transmembranePotentials;
  transmembranePotentials.reserve(nCompartments_);
  for (int i = 0; i < nCompartments_; i++)
  {
    transmembranePotentials.push_back(transmembranePotential_[i]);
  }

  std::vector<std::shared_ptr<FieldVariableType>> compartmentRelativeFactors;
  compartmentRelativeFactors.reserve(nCompartments_);
  for (int i = 0; i < nCompartments_; i++)
  {
    compartmentRelativeFactors.push_back(compartmentRelativeFactor_[i]);
  }

  return std::make_tuple(geometryField, this->fiberDirection_, this->flowPotential_, extraCellularPotential_,
                         transmembranePotentials, compartmentRelativeFactors, relativeFactorTotal_, activeStressTotal_);
}

} // namespace
