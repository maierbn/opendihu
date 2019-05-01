#include "data_management/specialized_solver/quasi_static_linear_elasticity.h"

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

template<typename DataLinearElasticityType>
QuasiStaticLinearElasticity<DataLinearElasticityType>::
QuasiStaticLinearElasticity(DihuContext context) :
  Data<FunctionSpaceType>::Data(context)
{
}

template<typename DataLinearElasticityType>
void QuasiStaticLinearElasticity<DataLinearElasticityType>::
initialize()
{
  // call initialize of base class
  Data<FunctionSpaceType>::initialize();
}

template<typename DataLinearElasticityType>
void QuasiStaticLinearElasticity<DataLinearElasticityType>::
createPetscObjects()
{
  LOG(DEBUG) << "QuasiStaticLinearElasticity::createPetscObject";

  assert(this->functionSpace_);

  this->activation_ = this->functionSpace_->template createFieldVariable<1>("activation");
  this->activeStress_ = this->functionSpace_->template createFieldVariable<9>("activeStress");
  this->strain_ = this->functionSpace_->template createFieldVariable<9>("strain");
  this->flowPotential_ = this->functionSpace_->template createFieldVariable<1>("flowPotential");
  this->fiberDirection_ = this->functionSpace_->template createFieldVariable<3>("fiberDirection");
}

template<typename DataLinearElasticityType>
std::shared_ptr<typename QuasiStaticLinearElasticity<DataLinearElasticityType>::FieldVariableType>
QuasiStaticLinearElasticity<DataLinearElasticityType>::
activation()
{
  return this->activation_;
}

template<typename DataLinearElasticityType>
std::shared_ptr<typename QuasiStaticLinearElasticity<DataLinearElasticityType>::StressFieldVariableType>
QuasiStaticLinearElasticity<DataLinearElasticityType>::
activeStress()
{
  return this->activeStress_;
}

template<typename DataLinearElasticityType>
std::shared_ptr<typename QuasiStaticLinearElasticity<DataLinearElasticityType>::StressFieldVariableType>
QuasiStaticLinearElasticity<DataLinearElasticityType>::
strain()
{
  return this->strain_;
}

template<typename DataLinearElasticityType>
std::shared_ptr<typename QuasiStaticLinearElasticity<DataLinearElasticityType>::GradientFieldVariableType>
QuasiStaticLinearElasticity<DataLinearElasticityType>::
fiberDirection()
{
  return this->fiberDirection_;
}

template<typename DataLinearElasticityType>
std::shared_ptr<typename QuasiStaticLinearElasticity<DataLinearElasticityType>::FieldVariableType>
QuasiStaticLinearElasticity<DataLinearElasticityType>::
flowPotential()
{
  return this->flowPotential_;
}

template<typename DataLinearElasticityType>
void QuasiStaticLinearElasticity<DataLinearElasticityType>::
setData(std::shared_ptr<DataLinearElasticityType> dataLinearElasticity)
{
  dataLinearElasticity_ = dataLinearElasticity;
}

template<typename DataLinearElasticityType>
void QuasiStaticLinearElasticity<DataLinearElasticityType>::
print()
{
  if (!VLOG_IS_ON(4))
    return;

  VLOG(4) << *this->activation_;
}

template<typename DataLinearElasticityType>
typename QuasiStaticLinearElasticity<DataLinearElasticityType>::OutputFieldVariables QuasiStaticLinearElasticity<DataLinearElasticityType>::
getOutputFieldVariables()
{
  return std::tuple_cat(
    dataLinearElasticity_->getOutputFieldVariables(),
    std::tuple<std::shared_ptr<FieldVariableType>>(this->activation_),
    std::tuple<std::shared_ptr<StressFieldVariableType>>(this->activeStress_),
    std::tuple<std::shared_ptr<StressFieldVariableType>>(this->strain_),
    std::tuple<std::shared_ptr<GradientFieldVariableType>>(this->fiberDirection_),
    std::tuple<std::shared_ptr<FieldVariableType>>(this->flowPotential_)
  );
}

} // namespace
