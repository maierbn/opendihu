#include "data_management/specialized_solver/reaction_diffusion_accelerator.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <numeric>
#include <memory>

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "control/dihu_context.h"
#include "utility/petsc_utility.h"

namespace Data
{

template<typename FunctionSpaceType,int nStates>
ReactionDiffusionAccelerator<FunctionSpaceType>::
ReactionDiffusionAccelerator(DihuContext context) :
  Data<FunctionSpaceType>::Data(context)
{
}

template<typename FunctionSpaceType,int nStates>
void ReactionDiffusionAccelerator<FunctionSpaceType>::
initialize()
{
  // call initialize of base class
  Data<FunctionSpaceType>::initialize();
}

template<typename FunctionSpaceType,int nStates>
void ReactionDiffusionAccelerator<FunctionSpaceType>::
createPetscObjects()
{
  LOG(DEBUG) << "ReactionDiffusionAccelerator::createPetscObject";

  assert(this->functionSpace_);

  this->states_ = this->functionSpace_->template createFieldVariable<nStates>("states");
  this->rates_ = this->functionSpace_->template createFieldVariable<nStates>("rates");
  this->intermediates_ = this->functionSpace_->template createFieldVariable<1>("intermediates");
  this->parameters_ = this->functionSpace_->template createFieldVariable<2>("parameters");
}

template<typename FunctionSpaceType,int nStates>
std::shared_ptr<typename ReactionDiffusionAccelerator<FunctionSpaceType>::FieldVariableType>
ReactionDiffusionAccelerator<FunctionSpaceType>::
states()
{
  return this->states_;
}

template<typename FunctionSpaceType,int nStates>
std::shared_ptr<typename ReactionDiffusionAccelerator<FunctionSpaceType>::StressFieldVariableType>
ReactionDiffusionAccelerator<FunctionSpaceType>::
rates()
{
  return this->rates_;
}

template<typename FunctionSpaceType,int nStates>
std::shared_ptr<typename ReactionDiffusionAccelerator<FunctionSpaceType>::StressFieldVariableType>
ReactionDiffusionAccelerator<FunctionSpaceType>::
intermediates()
{
  return this->intermediates_;
}

template<typename FunctionSpaceType,int nStates>
std::shared_ptr<typename ReactionDiffusionAccelerator<FunctionSpaceType>::VectorFieldVariableType>
ReactionDiffusionAccelerator<FunctionSpaceType>::
parameters()
{
  return this->parameters_;
}

template<typename FunctionSpaceType,int nStates>
typename ReactionDiffusionAccelerator<FunctionSpaceType>::OutputFieldVariables ReactionDiffusionAccelerator<FunctionSpaceType>::
getOutputFieldVariables()
{
  return std::make_tuple(this->functionSpace_->geometryField(), this->states_, this->rates_, this->intermediates_, this->parameters_);
}

} // namespace
