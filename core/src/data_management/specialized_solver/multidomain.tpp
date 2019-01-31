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
  nCompartments_ = nCompartments;

  // call initialize of base class
  Data<FunctionSpaceType>::initialize();
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

  assert(this->functionSpace_);

  for (int k = 0; k < nCompartments_; k++)
  {
    std::stringstream transmembranePotentialSolutionName;
    transmembranePotentialSolutionName << "Vm_solution_" << k;
    this->transmembranePotentialSolution_.push_back(this->functionSpace_->template createFieldVariable<1>(transmembranePotentialSolutionName.str()));

    std::stringstream transmembranePotentialName;
    transmembranePotentialName << "Vm_" << k;
    this->transmembranePotential_.push_back(this->functionSpace_->template createFieldVariable<1>(transmembranePotentialName.str()));

    std::stringstream compartmentRelativeFactorName;
    compartmentRelativeFactorName << "fr_" << k;
    this->compartmentRelativeFactor_.push_back(this->functionSpace_->template createFieldVariable<1>(compartmentRelativeFactorName.str()));
  }

  this->flowPotential_ = this->functionSpace_->template createFieldVariable<1>("flowPotential");
  this->fiberDirection_ = this->functionSpace_->template createFieldVariable<3>("fiberDirection");
  this->extraCellularPotential_ = this->functionSpace_->template createFieldVariable<1>("phi_e");
  this->zero_ = this->functionSpace_->template createFieldVariable<1>("zero");
  this->relativeFactorTotal_ = this->functionSpace_->template createFieldVariable<1>("relativeFactorTotal");
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
relativeFactorTotal()
{
  return this->relativeFactorTotal_;
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> Multidomain<FunctionSpaceType>::
zero()
{
  return this->zero_;
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
typename Multidomain<FunctionSpaceType>::OutputFieldVariables Multidomain<FunctionSpaceType>::
getOutputFieldVariables()
{
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
                         transmembranePotentials, compartmentRelativeFactors);
}

} // namespace
