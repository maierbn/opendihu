#include "data_management/specialized_solver/dynamic_hyperelasticity_solver.h"

namespace Data
{

template<typename FunctionSpaceType>
  DynamicHyperelasticitySolver<FunctionSpaceType>::
DynamicHyperelasticitySolver(DihuContext context) :
  Data<FunctionSpaceType>::Data(context)
{
}

template<typename FunctionSpaceType>
void DynamicHyperelasticitySolver<FunctionSpaceType>::
initialize(std::shared_ptr<DisplacementsFieldVariableType> displacements)
{
  // call initialize of base class, this calls createPetscObjects
  Data<FunctionSpaceType>::initialize();

  displacements_ = displacements;
}

template<typename FunctionSpaceType>
void DynamicHyperelasticitySolver<FunctionSpaceType>::
createPetscObjects()
{
  LOG(DEBUG) << "DynamicHyperelasticitySolver::createPetscObject";

  assert(this->functionSpace_);

  std::vector<std::string> displacementsComponentNames({"x","y","z"});
  velocity_                         = this->functionSpace_->template createFieldVariable<3>("v", displacementsComponentNames);
  acceleration_                     = this->functionSpace_->template createFieldVariable<3>("a", displacementsComponentNames);
  internalVirtualWorkCompressible_  = this->functionSpace_->template createFieldVariable<3>("∂Wint_compressible", displacementsComponentNames);
  displacementsCompressible_        = this->functionSpace_->template createFieldVariable<3>("u_compressible", displacementsComponentNames);
  velocityCompressible_             = this->functionSpace_->template createFieldVariable<3>("v_compressible", displacementsComponentNames);
  externalVirtualWork_              = this->functionSpace_->template createFieldVariable<3>("∂Wext", displacementsComponentNames);
}

//! field variable of u
template<typename FunctionSpaceType>
std::shared_ptr<typename DynamicHyperelasticitySolver<FunctionSpaceType>::DisplacementsFieldVariableType> DynamicHyperelasticitySolver<FunctionSpaceType>::
displacements()
{
  return this->displacements_;
}

//! field variable of v
template<typename FunctionSpaceType>
std::shared_ptr<typename DynamicHyperelasticitySolver<FunctionSpaceType>::DisplacementsFieldVariableType> DynamicHyperelasticitySolver<FunctionSpaceType>::
velocity()
{
  return this->velocity_;
}

//! field variable of a
template<typename FunctionSpaceType>
std::shared_ptr<typename DynamicHyperelasticitySolver<FunctionSpaceType>::DisplacementsFieldVariableType> DynamicHyperelasticitySolver<FunctionSpaceType>::
acceleration()
{
  return this->acceleration_;
}


//! field variable of u_compressible
template<typename FunctionSpaceType>
std::shared_ptr<typename DynamicHyperelasticitySolver<FunctionSpaceType>::DisplacementsFieldVariableType> DynamicHyperelasticitySolver<FunctionSpaceType>::
displacementsCompressible()
{
  return this->displacementsCompressible_;
}

//! field variable of v_compressible
template<typename FunctionSpaceType>
std::shared_ptr<typename DynamicHyperelasticitySolver<FunctionSpaceType>::DisplacementsFieldVariableType> DynamicHyperelasticitySolver<FunctionSpaceType>::
velocityCompressible()
{
  return this->velocityCompressible_;
}

//! field variable of ∂W_int_compressible
template<typename FunctionSpaceType>
std::shared_ptr<typename DynamicHyperelasticitySolver<FunctionSpaceType>::DisplacementsFieldVariableType> DynamicHyperelasticitySolver<FunctionSpaceType>::
internalVirtualWorkCompressible()
{
  return this->internalVirtualWorkCompressible_;
}

//! field variable of ∂Wext
template<typename FunctionSpaceType>
std::shared_ptr<typename DynamicHyperelasticitySolver<FunctionSpaceType>::DisplacementsFieldVariableType> DynamicHyperelasticitySolver<FunctionSpaceType>::
externalVirtualWork()
{
  return this->externalVirtualWork_;
}

template<typename FunctionSpaceType>
typename DynamicHyperelasticitySolver<FunctionSpaceType>::FieldVariablesForOutputWriter
DynamicHyperelasticitySolver<FunctionSpaceType>::
getFieldVariablesForOutputWriter()
{
  return std::make_tuple(
    std::shared_ptr<DisplacementsFieldVariableType>(std::make_shared<typename FunctionSpaceType::GeometryFieldType>(this->functionSpace_->geometryField())), // geometry
    std::shared_ptr<DisplacementsFieldVariableType>(this->displacements_),              // displacements_
    std::shared_ptr<DisplacementsFieldVariableType>(this->velocity_),
    std::shared_ptr<DisplacementsFieldVariableType>(this->acceleration_),
    std::shared_ptr<DisplacementsFieldVariableType>(this->displacementsCompressible_),
    std::shared_ptr<DisplacementsFieldVariableType>(this->velocityCompressible_),
    std::shared_ptr<DisplacementsFieldVariableType>(this->internalVirtualWorkCompressible_),
    std::shared_ptr<DisplacementsFieldVariableType>(this->externalVirtualWork_)
  );
}

} // namespace
