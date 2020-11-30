#include "data_management/specialized_solver/dummy.h"

namespace Data
{

template<typename FunctionSpaceType>
Dummy<FunctionSpaceType>::
Dummy(DihuContext context) :
  Data<FunctionSpaceType>::Data(context)
{
}

template<typename FunctionSpaceType>
void Dummy<FunctionSpaceType>::
initialize()
{
  // call initialize of base class, this calls createPetscObjects()
  Data<FunctionSpaceType>::initialize();

  // create the slot connector data object
  slotConnectorData_ = std::make_shared<SlotConnectorDataType>();

  // no field variables will be added
}

template<typename FunctionSpaceType>
void Dummy<FunctionSpaceType>::
createPetscObjects()
{
  // no field variables are created
}

template<typename FunctionSpaceType>
std::shared_ptr<typename Dummy<FunctionSpaceType>::SlotConnectorDataType> Dummy<FunctionSpaceType>::
getSlotConnectorData()
{
  // return the slot connector data object
  return this->slotConnectorData_;
}

} // namespace
