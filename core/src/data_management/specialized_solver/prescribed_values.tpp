#include "data_management/specialized_solver/prescribed_values.h"

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

template<typename FunctionSpaceType, int nComponents1, int nComponents2>
PrescribedValues<FunctionSpaceType,nComponents1,nComponents2>::
PrescribedValues(DihuContext context) :
  Data<FunctionSpaceType>::Data(context)
{
}

template<typename FunctionSpaceType, int nComponents1, int nComponents2>
void PrescribedValues<FunctionSpaceType,nComponents1,nComponents2>::
initialize(std::vector<std::string> fieldVariable1Names, std::vector<std::string> fieldVariable2Names)
{
  fieldVariable1Names_ = fieldVariable1Names;
  fieldVariable2Names_ = fieldVariable2Names;

  // call initialize of base class, this calls createPetscObjects()
  Data<FunctionSpaceType>::initialize();

  // create th output connector data object
  outputConnectorData_ = std::make_shared<OutputConnectorDataType>();

  // add all needed field variables to be transferred

  // add field variables with `nComponents1` components
  for (int fieldVariable1No = 0; fieldVariable1No < this->fieldVariables1_.size(); fieldVariable1No++)
  {
    // add component 0 of the field variable, the component 0 is hard-coded for now
    outputConnectorData_->addFieldVariable(this->fieldVariables1_[fieldVariable1No], 0);
  }

  // add field variables with `nComponents2` components
  for (int fieldVariable2No = 0; fieldVariable2No < this->fieldVariables2_.size(); fieldVariable2No++)
  {
    // add component 0 of the field variable, the component 0 is hard-coded for now
    outputConnectorData_->addFieldVariable2(this->fieldVariables2_[fieldVariable2No], 0);
  }

  LOG(DEBUG) << "create outputConnectorData: " << *this->outputConnectorData_;

  // parse slot names for all output connector data slots
  this->context_.getPythonConfig().getOptionVector("slotNames", outputConnectorData_->slotNames);
}

template<typename FunctionSpaceType, int nComponents1, int nComponents2>
void PrescribedValues<FunctionSpaceType,nComponents1,nComponents2>::
createPetscObjects()
{
  assert(this->functionSpace_);

  // create field variables with `nComponents1` components, using the given names
  for (int fieldVariable1No = 0; fieldVariable1No < fieldVariable1Names_.size(); fieldVariable1No++)
  {
    std::string name = fieldVariable1Names_[fieldVariable1No];
    this->fieldVariables1_.push_back(this->functionSpace_->template createFieldVariable<nComponents1>(name));
  }

  // create field variables with `nComponents2` components, using the given names
  for (int fieldVariable2No = 0; fieldVariable2No < fieldVariable2Names_.size(); fieldVariable2No++)
  {
    std::string name = fieldVariable2Names_[fieldVariable2No];
    this->fieldVariables2_.push_back(this->functionSpace_->template createFieldVariable<nComponents2>(name));
  }
}

template<typename FunctionSpaceType, int nComponents1, int nComponents2>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents1>> PrescribedValues<FunctionSpaceType,nComponents1,nComponents2>::
fieldVariable1(int index)
{
  assert(index < this->outputConnectorData_->variable1.size());
  return this->outputConnectorData_->variable1[index].values;
}

template<typename FunctionSpaceType, int nComponents1, int nComponents2>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents2>> PrescribedValues<FunctionSpaceType,nComponents1,nComponents2>::
fieldVariable2(int index)
{
  assert(index < this->outputConnectorData_->variable2.size());
  return this->outputConnectorData_->variable2[index].values;
}

template<typename FunctionSpaceType, int nComponents1, int nComponents2>
std::shared_ptr<typename PrescribedValues<FunctionSpaceType,nComponents1,nComponents2>::OutputConnectorDataType> PrescribedValues<FunctionSpaceType,nComponents1,nComponents2>::
getOutputConnectorData()
{
  // return the output connector data object
  return this->outputConnectorData_;
}

template<typename FunctionSpaceType, int nComponents1, int nComponents2>
typename PrescribedValues<FunctionSpaceType,nComponents1,nComponents2>::FieldVariablesForOutputWriter PrescribedValues<FunctionSpaceType,nComponents1,nComponents2>::
getFieldVariablesForOutputWriter()
{
  // these field variables will be written to output files by the output writer

  // get the geometry field, which is always needed, from the function space
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> geometryField
    = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,3>>(this->functionSpace_->geometryField());

  // update the pointer of the field variables, because the most recent field variable may be in outputConnectorData_
  // (it could have been changed during transfer)
  for (int fieldVariable1No = 0; fieldVariable1No < fieldVariables1_.size(); fieldVariable1No++)
  {
    this->fieldVariables1_[fieldVariable1No] = fieldVariable1(fieldVariable1No);
  }
  for (int fieldVariable2No = 0; fieldVariable2No < fieldVariables2_.size(); fieldVariable2No++)
  {
    this->fieldVariables2_[fieldVariable2No] = fieldVariable2(fieldVariable2No);
  }


  return std::make_tuple(
    geometryField,
    this->fieldVariables1_,
    this->fieldVariables2_
  );
}

} // namespace
