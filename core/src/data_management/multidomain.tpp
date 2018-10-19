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

template<typename FunctionSpaceType, int nStatesCellML>
Multidomain<FunctionSpaceType,nStatesCellML>::
Multidomain(DihuContext context) :
  Data<FunctionSpaceType>::Data(context)
{
}

template<typename FunctionSpaceType,int nStatesCellML>
void Multidomain<FunctionSpaceType,nStatesCellML>::
initialize(int nCompartments)
{
  nCompartments_ = nCompartments;

  // call initialize of base class
  Data<FunctionSpaceType>::initialize();
}

template<typename FunctionSpaceType,int nStatesCellML>
void Multidomain<FunctionSpaceType,nStatesCellML>::
createPetscObjects()
{
  LOG(DEBUG) << "Multidomain::createPetscObject for " << nCompartments_ << " compartments.";

  // create field variables that have one for every compartment
  subcellularStates_.reserve(nCompartments_);
  subcellularIncrement_.reserve(nCompartments_);
  ionicCurrent_.reserve(nCompartments_);
  transmembranePotential_.reserve(nCompartments_);
  compartmentRelativeFactor_.reserve(nCompartments_);

  assert(this->functionSpace_);

  for (int k = 0; k < nCompartments_; k++)
  {
    std::stringstream subcellularStatesName;
    subcellularStatesName << "y_" << k;
    this->subcellularStates_.push_back(this->functionSpace_->template createFieldVariable<nStatesCellML>(subcellularStatesName.str()));

    std::stringstream subcellularIncrementName;
    subcellularIncrementName << "Δy_" << k;
    this->subcellularIncrement_.push_back(this->functionSpace_->template createFieldVariable<nStatesCellML>(subcellularIncrementName.str()));

    std::stringstream ionicCurrentName;
    ionicCurrentName << "I_" << k;
    this->ionicCurrent_.push_back(this->functionSpace_->template createFieldVariable<1>(ionicCurrentName.str()));

    std::stringstream transmembranePotentialName;
    transmembranePotentialName << "Vm_" << k;
    this->transmembranePotential_.push_back(this->functionSpace_->template createFieldVariable<1>(transmembranePotentialName.str()));

    std::stringstream compartmentRelativeFactorName;
    compartmentRelativeFactorName << "fr_" << k;
    this->compartmentRelativeFactor_.push_back(this->functionSpace_->template createFieldVariable<1>(compartmentRelativeFactorName.str()));
  }

  this->flowPotential_ = this->functionSpace_->template createFieldVariable<1>("flowPotential");
  this->fibreDirection_ = this->functionSpace_->template createFieldVariable<3>("fibreDirection");
  this->extraCellularPotential_ = this->functionSpace_->template createFieldVariable<1>("phi_e");
  this->zero_ = this->functionSpace_->template createFieldVariable<1>("zero");
}

template<typename FunctionSpaceType,int nStatesCellML>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> Multidomain<FunctionSpaceType,nStatesCellML>::
fibreDirection()
{
  return this->fibreDirection_;
}

template<typename FunctionSpaceType,int nStatesCellML>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> Multidomain<FunctionSpaceType,nStatesCellML>::
flowPotential()
{
  return this->flowPotential_;
}

template<typename FunctionSpaceType,int nStatesCellML>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> Multidomain<FunctionSpaceType,nStatesCellML>::
extraCellularPotential()
{
  return this->extraCellularPotential_;
}

template<typename FunctionSpaceType,int nStatesCellML>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> Multidomain<FunctionSpaceType,nStatesCellML>::
transmembranePotential(int compartmentNo)
{
  assert(compartmentNo >= 0 && compartmentNo < nCompartments_);
  return this->transmembranePotential_[compartmentNo];
}

template<typename FunctionSpaceType,int nStatesCellML>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> Multidomain<FunctionSpaceType,nStatesCellML>::
compartmentRelativeFactor(int compartmentNo)
{
  assert(compartmentNo >= 0 && compartmentNo < nCompartments_);
  return this->compartmentRelativeFactor_[compartmentNo];
}

template<typename FunctionSpaceType,int nStatesCellML>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nStatesCellML>> Multidomain<FunctionSpaceType,nStatesCellML>::
subcellularStates(int compartmentNo)
{
  assert(compartmentNo >= 0 && compartmentNo < nCompartments_);
  return this->subcellularStates_[compartmentNo];
}

template<typename FunctionSpaceType,int nStatesCellML>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nStatesCellML>> Multidomain<FunctionSpaceType,nStatesCellML>::
subcellularIncrement(int compartmentNo)
{
  assert(compartmentNo >= 0 && compartmentNo < nCompartments_);
  return this->subcellularIncrement_[compartmentNo];
}

template<typename FunctionSpaceType,int nStatesCellML>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> Multidomain<FunctionSpaceType,nStatesCellML>::
ionicCurrent(int compartmentNo)
{
  assert(compartmentNo >= 0 && compartmentNo < nCompartments_);
  return this->ionicCurrent_[compartmentNo];
}

template<typename FunctionSpaceType,int nStatesCellML>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> Multidomain<FunctionSpaceType,nStatesCellML>::
zero()
{
  return this->zero_;
}

template<typename FunctionSpaceType,int nStatesCellML>
void Multidomain<FunctionSpaceType,nStatesCellML>::
print() // use override in stead of extending the parents' print output.This way "solution" is still in the end.
{
  if (!VLOG_IS_ON(4))
    return;

  VLOG(4) << *this->fibreDirection_;
}

template<typename FunctionSpaceType,int nStatesCellML>
typename Multidomain<FunctionSpaceType,nStatesCellML>::OutputFieldVariables Multidomain<FunctionSpaceType,nStatesCellML>::
getOutputFieldVariables()
{
  std::vector<std::shared_ptr<FieldVariableType>> transmembranePotentials;
  transmembranePotentials.reserve(nCompartments_);
  for (int i = 0; i < nCompartments_; i++)
  {
    transmembranePotentials.push_back(transmembranePotential_[i]);
  }

  std::vector<std::shared_ptr<CellMLFieldVariableType>> subcellularStates;
  subcellularStates.reserve(nCompartments_);
  for (int i = 0; i < nCompartments_; i++)
  {
    subcellularStates.push_back(subcellularStates_[i]);
  }

  std::vector<std::shared_ptr<FieldVariableType>> compartmentRelativeFactors;
  compartmentRelativeFactors.reserve(nCompartments_);
  for (int i = 0; i < nCompartments_; i++)
  {
    compartmentRelativeFactors.push_back(compartmentRelativeFactor_[i]);
  }

  return std::make_tuple(this->fibreDirection_, extraCellularPotential_, transmembranePotentials, subcellularStates, compartmentRelativeFactors);
}

} // namespace
