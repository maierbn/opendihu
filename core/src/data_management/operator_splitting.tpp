#include "data_management/operator_splitting.h"

namespace Data
{

template<typename TimeStepping1, typename TimeStepping2>
OperatorSplitting<TimeStepping1,TimeStepping2>::
OperatorSplitting(DihuContext context) :
  Data<typename TimeStepping1::FunctionSpace>(context)
{
}

template<typename TimeStepping1, typename TimeStepping2>
void OperatorSplitting<TimeStepping1,TimeStepping2>::
initialize(TimeStepping1 &timeStepping1, TimeStepping2 &timeStepping2)
{
  Data<typename TimeStepping1::FunctionSpace>::initialize();

  // create the outputConnectorData_ object and assign the two objects of the timestepping schemes
  outputConnectorData_ = std::make_shared<OutputConnectorDataType>();

  std::get<0>(*outputConnectorData_) = timeStepping1.getOutputConnectorData();
  std::get<1>(*outputConnectorData_) = timeStepping2.getOutputConnectorData();
}

template<typename TimeStepping1, typename TimeStepping2>
std::shared_ptr<typename OperatorSplitting<TimeStepping1,TimeStepping2>::OutputConnectorDataType>
OperatorSplitting<TimeStepping1,TimeStepping2>::
getOutputConnectorData()
{
  return outputConnectorData_;
}

} // namespace
