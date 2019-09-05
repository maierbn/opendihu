#include "data_management/cellml_adapter.h"

namespace Data
{

template <int nStates, int nIntermediates, typename FunctionSpaceType>
CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::
CellmlAdapter(DihuContext context) :
  Data<FunctionSpaceType>::Data(context)
{
}

template <int nStates, int nIntermediates, typename FunctionSpaceType>
void CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::
initialize()
{
  // call initialize of base class
  Data<FunctionSpaceType>::initialize();
}

template <int nStates, int nIntermediates, typename FunctionSpaceType>
void CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::
setStatesVariable(std::shared_ptr<CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::FieldVariableStates> states)
{
  this->states_ = states;
}

template <int nStates, int nIntermediates, typename FunctionSpaceType>
void CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::
createPetscObjects()
{
  this->intermediates_ = this->functionSpace_->template createFieldVariable<nIntermediates>("intermediates");
  this->intermediates_->setRepresentationContiguous();
}

template <int nStates, int nIntermediates, typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nIntermediates>> CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::
intermediates()
{
  return this->intermediates_;
}

template <int nStates, int nIntermediates, typename FunctionSpaceType>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nStates>> CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::
states()
{
  return this->states_;
}

template <int nStates, int nIntermediates, typename FunctionSpaceType>
typename CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::FieldVariablesForOutputWriter CellmlAdapter<nStates,nIntermediates,FunctionSpaceType>::
getFieldVariablesForOutputWriter()
{
  // these field variables will be written to output files
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> geometryField
    = std::make_shared<FieldVariable::FieldVariable<FunctionSpaceType,3>>(this->functionSpace_->geometryField());

  return std::make_tuple(geometryField, intermediates_, states_);
}

} // namespace
