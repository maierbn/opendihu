#include "time_stepping_scheme/time_stepping_scheme_ode_transferable_solution_data.h"

#include <Python.h>  // has to be the first included header
#include <vector>

#include "utility/python_utility.h"

namespace TimeSteppingScheme
{

template<typename FunctionSpaceType, int nComponents, typename DiscretizableInTimeType>
typename TimeSteppingSchemeOdeOutputConnectorDataType<FunctionSpaceType, nComponents, DiscretizableInTimeType>::OutputConnectorDataType
TimeSteppingSchemeOdeOutputConnectorDataType<FunctionSpaceType, nComponents, DiscretizableInTimeType>::
getOutputConnectorData()
{
  return this->data_->getOutputConnectorData();
}

//! output the given data for debugging
template<typename FunctionSpaceType, int nComponents, typename DiscretizableInTimeType>
std::string TimeSteppingSchemeOdeOutputConnectorDataType<FunctionSpaceType, nComponents, DiscretizableInTimeType>::
getString(typename TimeSteppingSchemeOdeOutputConnectorDataType<FunctionSpaceType, nComponents, DiscretizableInTimeType>::OutputConnectorDataType &data)
{
  return this->data_->getString(data);
}


} // namespace
