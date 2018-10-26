#include "operator_splitting/solution_vector_mapping/solution_vector_mapping_finite_element_method_cellml.h"

#include <vector>
#include <tuple>
#include "easylogging++.h"

template<typename FunctionSpaceType1, typename FunctionSpaceType2, int nComponents2>
void SolutionVectorMapping<
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>>,   // <fieldVariableType,componentNo>
  std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>>, int, double>
>::transfer(const std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>> &transferableSolutionData1,
            const std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>>, int, double> &transferableSolutionData2)
{
  // call the transfer function of the parent class
  SolutionVectorMapping<
    std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>>, int, double>,
    std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>>, int, double>
  >::transfer(std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>>, int, double>(transferableSolutionData1, 0, 1.0),
              transferableSolutionData2);
}

template<typename FunctionSpaceType1, typename FunctionSpaceType2, int nComponents2>
void SolutionVectorMapping<
  std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>>, int, double>,
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>>   // <fieldVariableType,componentNo>
>::transfer(const std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>>, int, double> &transferableSolutionData2,
            const std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>> &transferableSolutionData1)
{
  // call the transfer function of the parent class
  SolutionVectorMapping<
    std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType2,nComponents2>>, int, double>,
    std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>>, int, double>
  >::transfer(transferableSolutionData2,
              std::tuple<std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType1,1>>, int, double>(transferableSolutionData1, 0, 1.0));
}
