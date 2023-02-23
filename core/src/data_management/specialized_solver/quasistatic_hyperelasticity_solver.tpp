#include "data_management/specialized_solver/quasistatic_hyperelasticity_solver.h"

namespace Data
{

template<typename FunctionSpaceType>
  QuasistaticHyperelasticitySolver<FunctionSpaceType>::
QuasistaticHyperelasticitySolver(DihuContext context) :
  Data<FunctionSpaceType>::Data(context)
{
}

template<typename FunctionSpaceType>
void QuasistaticHyperelasticitySolver<FunctionSpaceType>::
createPetscObjects()
{
  LOG(DEBUG) << "QuasistaticHyperelasticitySolver::createPetscObjects";

  assert(this->functionSpace_);

  std::vector<std::string> displacementsComponentNames({"x","y","z"});
  displacements_                    = this->functionSpace_->template createFieldVariable<3>("u", displacementsComponentNames);
  velocities_                       = this->functionSpace_->template createFieldVariable<3>("v", displacementsComponentNames);
  internalVirtualWork_              = this->functionSpace_->template createFieldVariable<3>("δWint_displacement", displacementsComponentNames);
  accelerationTerm_                 = this->functionSpace_->template createFieldVariable<3>("δWint_acceleration", displacementsComponentNames);
  externalVirtualWorkDead_          = this->functionSpace_->template createFieldVariable<3>("δWext", displacementsComponentNames);
}

//! field variable of u
template<typename FunctionSpaceType>
std::shared_ptr<typename QuasistaticHyperelasticitySolver<FunctionSpaceType>::DisplacementsFieldVariableType> QuasistaticHyperelasticitySolver<FunctionSpaceType>::
displacements()
{
  return this->displacements_;
}

//! field variable of v
template<typename FunctionSpaceType>
std::shared_ptr<typename QuasistaticHyperelasticitySolver<FunctionSpaceType>::DisplacementsFieldVariableType> QuasistaticHyperelasticitySolver<FunctionSpaceType>::
velocities()
{
  return this->velocities_;
}

//! field variable of u_compressible
template<typename FunctionSpaceType>
std::shared_ptr<typename QuasistaticHyperelasticitySolver<FunctionSpaceType>::DisplacementsFieldVariableType> QuasistaticHyperelasticitySolver<FunctionSpaceType>::
externalVirtualWorkDead()
{
  return this->externalVirtualWorkDead_;
}

//! field variable of ∂W_int_compressible
template<typename FunctionSpaceType>
std::shared_ptr<typename QuasistaticHyperelasticitySolver<FunctionSpaceType>::DisplacementsFieldVariableType> QuasistaticHyperelasticitySolver<FunctionSpaceType>::
internalVirtualWork()
{
  return this->internalVirtualWork_;
}

//! field variable of
template<typename FunctionSpaceType>
std::shared_ptr<typename QuasistaticHyperelasticitySolver<FunctionSpaceType>::DisplacementsFieldVariableType> QuasistaticHyperelasticitySolver<FunctionSpaceType>::
accelerationTerm()
{
  return this->accelerationTerm_;
}

template<typename FunctionSpaceType>
typename QuasistaticHyperelasticitySolver<FunctionSpaceType>::FieldVariablesForOutputWriter
QuasistaticHyperelasticitySolver<FunctionSpaceType>::
getFieldVariablesForOutputWriter()
{
  return std::make_tuple(
    std::shared_ptr<DisplacementsFieldVariableType>(std::make_shared<typename FunctionSpaceType::GeometryFieldType>(this->functionSpace_->geometryField())), // geometry
    std::shared_ptr<DisplacementsFieldVariableType>(this->displacements_),              // displacements_
    std::shared_ptr<DisplacementsFieldVariableType>(this->velocities_),
    std::shared_ptr<DisplacementsFieldVariableType>(this->internalVirtualWork_),
    std::shared_ptr<DisplacementsFieldVariableType>(this->accelerationTerm_),
    std::shared_ptr<DisplacementsFieldVariableType>(this->externalVirtualWorkDead_)
  );
}

} // namespace
