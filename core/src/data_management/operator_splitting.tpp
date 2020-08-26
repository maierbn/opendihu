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

  // create the slotConnectorData_ object and assign the two objects of the timestepping schemes
  slotConnectorData_ = std::make_shared<SlotConnectorDataType>();

  std::get<0>(*slotConnectorData_) = timeStepping1.getSlotConnectorData();
  std::get<1>(*slotConnectorData_) = timeStepping2.getSlotConnectorData();

  timeStepping1_ = std::make_shared<TimeStepping1>(timeStepping1);
  timeStepping2_ = std::make_shared<TimeStepping2>(timeStepping2);
}

template<typename TimeStepping1, typename TimeStepping2>
std::shared_ptr<typename OperatorSplitting<TimeStepping1,TimeStepping2>::SlotConnectorDataType>
OperatorSplitting<TimeStepping1,TimeStepping2>::
getSlotConnectorData()
{
  // get the slot connector data again, such that prepareForGetSlotConnectorData() will be called
  std::get<0>(*slotConnectorData_) = timeStepping1_->getSlotConnectorData();
  std::get<1>(*slotConnectorData_) = timeStepping2_->getSlotConnectorData();

  return slotConnectorData_;
}

//! get pointers to all field variables that can be written by output writers
template<typename TimeStepping1, typename TimeStepping2>
typename OperatorSplitting<TimeStepping1,TimeStepping2>::FieldVariablesForOutputWriter OperatorSplitting<TimeStepping1,TimeStepping2>::
getFieldVariablesForOutputWriter()
{
  return timeStepping1_->data().getFieldVariablesForOutputWriter();
}


} // namespace
