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
  LOG(DEBUG) << "Multidomain::createPetscObject for " << nCompartments_ << " compartments.";
  transmembranePotential_.reserve(nCompartments);
  for (int i = 0; i < nCompartments; i++)
  {
    std::stringstream fieldVariableName;
    fieldVariableName << "Vm_" << i;
    transmembranePotential_.emplace_back(std::make_shared<FieldVariableType>(this->functionSpace_->template createFieldVariable<1>(fieldVariableName.str())));
  }

  this->fibreDirection_ = this->functionSpace_->template createFieldVariable<3>("fibreDirection");
  this->extraCellularPotential_ = this->functionSpace_->template createFieldVariable<1>("phi_e");
}

template<typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> Multidomain<FunctionSpaceType>::
fibreDirection()
{
  return this->fibreDirection_;
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
void Multidomain<FunctionSpaceType>::
print() // use override in stead of extending the parents' print output.This way "solution" is still in the end.
{
  if (!VLOG_IS_ON(4))
    return;

  VLOG(4) << *this->fibreDirection_;
}

template<typename FunctionSpaceType>
typename MultipleInstances<FunctionSpaceType,BaseDataType>::OutputFieldVariables Multidomain<FunctionSpaceType>::
getOutputFieldVariables()
{
  std::vector<FunctionSpaceType> outputFieldVariables;
  outputFieldVariables.reserve(nCompartments_);

  for (int i = 0; i < nCompartments_; i++)
  {
    outputFieldVariables.push_back(*transmembranePotential_[i]);
  }

  return std::make_tuple(this->fibreDirection_, extraCellularPotential_, outputFieldVariables);
}

} // namespace
