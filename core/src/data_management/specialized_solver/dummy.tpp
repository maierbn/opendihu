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

  // create the output connector data object
  outputConnectorData_ = std::make_shared<OutputConnectorDataType>();

  // no field variables will be added
}

template<typename FunctionSpaceType>
void Dummy<FunctionSpaceType>::
createPetscObjects()
{
  // no field variables are created
}

template<typename FunctionSpaceType>
std::shared_ptr<typename Dummy<FunctionSpaceType>::OutputConnectorDataType> Dummy<FunctionSpaceType>::
getOutputConnectorData()
{
  // return the output connector data object
  return this->outputConnectorData_;
}

} // namespace
