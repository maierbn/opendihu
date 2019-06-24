#include "time_stepping_scheme/time_stepping_scheme_ode_transferable_solution_data.h"

#include <Python.h>  // has to be the first included header
#include <vector>

#include "utility/python_utility.h"

namespace TimeSteppingScheme
{

template<typename FunctionSpaceType, int nComponents, typename DiscretizableInTimeType>
typename TimeSteppingSchemeOdeTransferableSolutionData<FunctionSpaceType, nComponents, DiscretizableInTimeType>::TransferableSolutionDataType
TimeSteppingSchemeOdeTransferableSolutionData<FunctionSpaceType, nComponents, DiscretizableInTimeType>::
getSolutionForTransfer()
{
  return this->data_->getSolutionForTransfer();
}

//! output the given data for debugging
template<typename FunctionSpaceType, int nComponents, typename DiscretizableInTimeType>
std::string TimeSteppingSchemeOdeTransferableSolutionData<FunctionSpaceType, nComponents, DiscretizableInTimeType>::
getString(typename TimeSteppingSchemeOdeTransferableSolutionData<FunctionSpaceType, nComponents, DiscretizableInTimeType>::TransferableSolutionDataType &data)
{
  return this->data_->getString(data);
}


} // namespace
