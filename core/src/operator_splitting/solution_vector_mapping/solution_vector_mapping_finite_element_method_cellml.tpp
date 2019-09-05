#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_finite_element_method_cellml.h"

#include <vector>
#include <tuple>
#include "easylogging++.h"

template<typename FunctionSpaceType1, typename OutputConnectorDataType2>
void SolutionVectorMapping<
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>>,   // <fieldVariableType,componentNo>
  OutputConnectorDataType2
>::transfer(const std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>> &transferableSolutionData1,
            const OutputConnectorDataType2 &transferableSolutionData2,
            const std::string transferSlotName)
{
  Data::ScaledFieldVariableComponent<FunctionSpaceType1,1> scaledFieldVariableComponent1;
  scaledFieldVariableComponent1.values = transferableSolutionData1;
  scaledFieldVariableComponent1.componentNo = 0;
  scaledFieldVariableComponent1.scalingFactor = 1.0;

  // call the transfer function of the parent class
  SolutionVectorMapping<
    Data::ScaledFieldVariableComponent<FunctionSpaceType1,1>,
    OutputConnectorDataType2
  >::transfer(scaledFieldVariableComponent1, transferableSolutionData2, transferSlotName);
}

template<typename OutputConnectorDataType1, typename FunctionSpaceType2>
void SolutionVectorMapping<
  OutputConnectorDataType1,
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,1>>   // <fieldVariableType,componentNo>
>::transfer(const OutputConnectorDataType1 &transferableSolutionData1,
            const std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,1>> &transferableSolutionData2,
            const std::string transferSlotName)
{
  Data::ScaledFieldVariableComponent<FunctionSpaceType2,1> scaledFieldVariableComponent2;
  scaledFieldVariableComponent2.values = transferableSolutionData2;
  scaledFieldVariableComponent2.componentNo = 0;
  scaledFieldVariableComponent2.scalingFactor = 1.0;

  // call the transfer function of the parent class
  SolutionVectorMapping<
    OutputConnectorDataType1,
    Data::ScaledFieldVariableComponent<FunctionSpaceType2,1>
  >::transfer(transferableSolutionData1, scaledFieldVariableComponent2, transferSlotName);
}

template<typename FunctionSpaceType1, typename FunctionSpaceType2>
void SolutionVectorMapping<
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>>,
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,1>>   // <fieldVariableType,componentNo>
>::transfer(const std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>> &transferableSolutionData1,
            const std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,1>> &transferableSolutionData2,
            const std::string transferSlotName)
{

  Data::ScaledFieldVariableComponent<FunctionSpaceType1,1> scaledFieldVariableComponent1;
  scaledFieldVariableComponent1.values = transferableSolutionData1;
  scaledFieldVariableComponent1.componentNo = 0;
  scaledFieldVariableComponent1.scalingFactor = 1.0;

  Data::ScaledFieldVariableComponent<FunctionSpaceType2,1> scaledFieldVariableComponent2;
  scaledFieldVariableComponent2.values = transferableSolutionData2;
  scaledFieldVariableComponent2.componentNo = 0;
  scaledFieldVariableComponent2.scalingFactor = 1.0;

  // call the transfer function of the parent class
  SolutionVectorMapping<
    Data::ScaledFieldVariableComponent<FunctionSpaceType1,1>,
    Data::ScaledFieldVariableComponent<FunctionSpaceType2,1>
  >::transfer(scaledFieldVariableComponent1, scaledFieldVariableComponent2, transferSlotName);
}
