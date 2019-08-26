#include "time_stepping_scheme/time_stepping_scheme_ode.h"

#include <Python.h>  // has to be the first included header
#include <vector>

#include "utility/python_utility.h"

namespace TimeSteppingScheme
{
// ----------------------
// specialization for CellML adapter
template<int nStates,int nIntermediates,typename FunctionSpaceType>
void TimeSteppingSchemeOde<CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>>::
initialize()
{
  TimeSteppingSchemeOdeBaseDiscretizable<CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>>::initialize();
  double prefactor = this->discretizableInTime_.prefactor();
  int outputComponentNo = this->discretizableInTime_.outputStateIndex();

  LOG(DEBUG) << "set CellML prefactor=" << prefactor << ", outputComponentNo=" << outputComponentNo;

  this->data_->setPrefactor(prefactor);
  this->data_->setOutputComponentNo(outputComponentNo);
}

template<int nStates,int nIntermediates,typename FunctionSpaceType>
typename TimeSteppingSchemeOde<CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>>::TransferableSolutionData
TimeSteppingSchemeOde<CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>>::
getSolutionForTransfer()
{
  // get the intermediates field variable with outputIntermediateIndex_
  std::tuple<std::shared_ptr<typename CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::FieldVariableTypeIntermediates>,int> intermediatesData(
    this->discretizableInTime_.intermediates(), this->discretizableInTime_.outputIntermediateIndex());

  // create pair of solution and intermediates data
  return std::make_pair<
    typename Data::TimeStepping<FunctionSpaceType,nStates>::TransferableSolutionDataType,
    std::tuple<std::shared_ptr<typename CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::FieldVariableTypeIntermediates>,int>
  >(
    this->data_->getSolutionForTransfer(),
    std::move(intermediatesData)
  );
}

//! output the given data for debugging
template<int nStates,int nIntermediates,typename FunctionSpaceType>
std::string TimeSteppingSchemeOde<CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>>::
getString(typename TimeSteppingSchemeOde<CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>>::TransferableSolutionData &data)
{
  return std::string("");
}

} // namespace
