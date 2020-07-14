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

  // initialize the output connector slots
  outputConnectorData_ = std::make_shared<OutputConnectorDataType>();

  // there is only one slot: the activation field variable
  outputConnectorData_->addFieldVariable(this->activation_);
  outputConnectorData_->addGeometryField(std::make_shared<VectorFieldVariableType>(this->functionSpace_->geometryField()));

  // parse slot names for all output connector data slots
  this->context_.getPythonConfig().getOptionVector("slotNames", outputConnectorData_->slotNames);
}

template<typename DataLinearElasticityType>
void QuasiStaticLinearElasticity<DataLinearElasticityType>::
createPetscObjects()
{
  LOG(DEBUG) << "QuasiStaticLinearElasticity::createPetscObject";

  assert(this->functionSpace_);

  std::vector<std::string> componentNames({"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"});

  // create all field variables that are needed
  this->activation_ = this->functionSpace_->template createFieldVariable<1>("activation");
  this->activeStress_ = this->functionSpace_->template createFieldVariable<9>("activeStress", componentNames);
  this->strain_ = this->functionSpace_->template createFieldVariable<9>("strain", componentNames);
  this->flowPotential_ = this->functionSpace_->template createFieldVariable<1>("flowPotential");
  this->rightHandSideActive_ = this->functionSpace_->template createFieldVariable<3>("rightHandSideActive");
  this->fiberDirection_ = this->functionSpace_->template createFieldVariable<3>("fiberDirection");
}

template<typename DataLinearElasticityType>
std::shared_ptr<typename QuasiStaticLinearElasticity<DataLinearElasticityType>::FieldVariableType>
QuasiStaticLinearElasticity<DataLinearElasticityType>::
activation()
{
  // get the most recent pointer to the field variable from the output connector slot, it may have changed because of sharing of the field variable
  this->activation_ = this->outputConnectorData_->variable1[0].values;
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
std::shared_ptr<typename QuasiStaticLinearElasticity<DataLinearElasticityType>::VectorFieldVariableType>
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
std::shared_ptr<typename QuasiStaticLinearElasticity<DataLinearElasticityType>::VectorFieldVariableType>
QuasiStaticLinearElasticity<DataLinearElasticityType>::
rightHandSideActive()
{
  return this->rightHandSideActive_;
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
void QuasiStaticLinearElasticity<DataLinearElasticityType>::
debug()
{
  return;
  std::shared_ptr<VectorFieldVariableType> solution = this->dataLinearElasticity_->solution();

  int nValues = solution->nDofsLocalWithoutGhosts();
  std::vector<Vec3> values;
  solution->getValuesWithoutGhosts(values);

  std::vector<Vec3> geometryValues;
  solution->functionSpace()->geometryField().getValuesWithoutGhosts(geometryValues);

  static int aa=0;

  for (int i = 0; i < nValues; i++)
  {
    std::stringstream s;
    s << values[i][0] << "," << values[i][1] << "," << values[i][2];
    values[i][2] = 1e-3;
    values[i][1] = 1e-3;
    values[i][0] = 1e-3;
    LOG(INFO) << "i: " << i << ", geometry: " << geometryValues[i] << ", solution: " << s.str() << " -> " << values[i];
  }
  aa++;
  LOG(INFO) << "--";
  solution->setValuesWithoutGhosts(values);
}


template<typename DataLinearElasticityType>
std::shared_ptr<typename QuasiStaticLinearElasticity<DataLinearElasticityType>::OutputConnectorDataType>
QuasiStaticLinearElasticity<DataLinearElasticityType>::
getOutputConnectorData()
{
  return outputConnectorData_;
}

template<typename DataLinearElasticityType>
typename QuasiStaticLinearElasticity<DataLinearElasticityType>::FieldVariablesForOutputWriter QuasiStaticLinearElasticity<DataLinearElasticityType>::
getFieldVariablesForOutputWriter()
{
  // these field variables will be written to output files
  return std::tuple_cat(
    dataLinearElasticity_->getFieldVariablesForOutputWriter(),
    std::tuple<std::shared_ptr<FieldVariableType>>(this->activation_),
    std::tuple<std::shared_ptr<StressFieldVariableType>>(this->activeStress_),
    std::tuple<std::shared_ptr<StressFieldVariableType>>(this->strain_),
    std::tuple<std::shared_ptr<VectorFieldVariableType>>(this->rightHandSideActive_),
    std::tuple<std::shared_ptr<VectorFieldVariableType>>(this->fiberDirection_),
    std::tuple<std::shared_ptr<FieldVariableType>>(this->flowPotential_)
  );
}

} // namespace
