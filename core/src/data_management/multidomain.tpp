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
  transmembranePotential_.reserve(nCompartments_);
  transmembranePotentialNextTimeStep_.reserve(nCompartments_);
  transmembraneIncrementNextTimeStep_.reserve(nCompartments_);
  ionicCurrentNextTimestep_.reserve(nCompartments_);

  assert(this->functionSpace_);

  for (int i = 0; i < nCompartments_; i++)
  {
    std::stringstream transmembranePotentialName;
    transmembranePotentialName << "Vm_" << i;
    this->transmembranePotential_.push_back(this->functionSpace_->template createFieldVariable<nStatesCellML>(transmembranePotentialName.str()));

    std::stringstream transmembranePotentialNextTimeStepName;
    transmembranePotentialNextTimeStepName << "Vm^{i+1}_" << i;
    this->transmembranePotentialNextTimeStep_.push_back(this->functionSpace_->template createFieldVariable<nStatesCellML>(transmembranePotentialNextTimeStepName.str()));

    std::stringstream transmembraneIncrementNextTimeStepName;
    transmembraneIncrementNextTimeStepName << "ΔVm^{i+1}_" << i;
    this->transmembraneIncrementNextTimeStep_.push_back(this->functionSpace_->template createFieldVariable<nStatesCellML>(transmembraneIncrementNextTimeStepName.str()));

    std::stringstream ionicCurrentNextTimestepName;
    ionicCurrentNextTimestepName << "I^{i+1}_" << i;
    this->ionicCurrentNextTimestep_.push_back(this->functionSpace_->template createFieldVariable<1>(ionicCurrentNextTimestepName.str()));
  }

  this->fibreDirection_ = this->functionSpace_->template createFieldVariable<3>("fibreDirection");
  this->extraCellularPotential_ = this->functionSpace_->template createFieldVariable<1>("phi_e");
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
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nStatesCellML>> Multidomain<FunctionSpaceType,nStatesCellML>::
transmembranePotential(int compartmentNo)
{
  assert(compartmentNo >= 0 && compartmentNo < nCompartments_);
  return this->transmembranePotential_[compartmentNo];
}

template<typename FunctionSpaceType,int nStatesCellML>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nStatesCellML>> Multidomain<FunctionSpaceType,nStatesCellML>::
transmembranePotentialNextTimeStep(int compartmentNo)
{
  assert(compartmentNo >= 0 && compartmentNo < nCompartments_);
  return this->transmembranePotentialNextTimeStep_[compartmentNo];
}

template<typename FunctionSpaceType,int nStatesCellML>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nStatesCellML>> Multidomain<FunctionSpaceType,nStatesCellML>::
transmembraneIncrementNextTimeStep(int compartmentNo)
{
  assert(compartmentNo >= 0 && compartmentNo < nCompartments_);
  return this->transmembraneIncrementNextTimeStep_[compartmentNo];
}

template<typename FunctionSpaceType,int nStatesCellML>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,1>> Multidomain<FunctionSpaceType,nStatesCellML>::
ionicCurrentNextTimeStep(int compartmentNo)
{
  assert(compartmentNo >= 0 && compartmentNo < nCompartments_);
  return this->ionicCurrentNextTimestep_[compartmentNo];
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
  std::vector<std::shared_ptr<CellMLFieldVariableType>> outputFieldVariables;
  outputFieldVariables.reserve(nCompartments_);

  for (int i = 0; i < nCompartments_; i++)
  {
    outputFieldVariables.push_back(transmembranePotential_[i]);
  }

  return std::make_tuple(this->fibreDirection_, extraCellularPotential_, outputFieldVariables);
}

} // namespace
