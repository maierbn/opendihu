#include "data_management/specialized_solver/quasi_static_hyperelasticity.h"

namespace Data
{

template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
  QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::
QuasiStaticHyperelasticity(DihuContext context) :
  Data<DisplacementsFunctionSpace>::Data(context)
{
}

template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
void QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::
initialize()
{
  // call initialize of base class
  Data<DisplacementsFunctionSpace>::initialize();
}

template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
void QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::
createPetscObjects()
{
  LOG(DEBUG) << "QuasiStaticHyperelasticity::createPetscObject";

  assert(this->displacementsFunctionSpace_);
  assert(this->pressureFunctionSpace_);
  assert(this->functionSpace_);


  std::vector<std::string> displacementsComponentNames({"x","y","z"});
  displacementsIncrement_ = this->displacementsFunctionSpace_->template createFieldVariable<3>("Δu", displacementsComponentNames);     //< Δu, the displacements increment that is solved for in the Newton scheme
  pressureIncrement_      = this->pressureFunctionSpace_->template createFieldVariable<1>("Δp");     //<  Δp, the pressure increment that is solved for in the NEwton scheme
  displacements_          = this->displacementsFunctionSpace_->template createFieldVariable<3>("u", displacementsComponentNames);     //< u, the displacements
  pressure_               = this->pressureFunctionSpace_->template createFieldVariable<1>("p");     //<  p, the pressure variable

  std::vector<std::string> componentNames({"S_11", "S_22", "S_33", "S_12", "S_13", "S_23"});
  pK2Stress_              = this->displacementsFunctionSpace_->template createFieldVariable<6>("PK2-Stress (Voigt)", componentNames);     //<  the symmetric PK2 stress tensor in Voigt notation
}


//! field variable of Δu
template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
std::shared_ptr<typename QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::DisplacementsFieldVariableType> QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::
displacementsIncrement()
{
  return this->displacementsIncrement_;
}

//! field variable of Δp
template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
std::shared_ptr<typename QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::PressureFieldVariableType> QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::
pressureIncrement()
{
  return this->pressureIncrement_;
}

//! field variable of u
template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
std::shared_ptr<typename QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::DisplacementsFieldVariableType> QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::
displacements()
{
  return this->displacements_;
}

//! field variable of p
template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
std::shared_ptr<typename QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::PressureFieldVariableType> QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::
pressure()
{
  return this->pressure_;
}

//! field variable of S
template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
std::shared_ptr<typename QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::StressFieldVariableType> QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::
pK2Stress()
{
  return this->pK2Stress_;
}


//! set the function space object that discretizes the pressure field variable
template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
void QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::
setPressureFunctionSpace(std::shared_ptr<PressureFunctionSpace> pressureFunctionSpace)
{
  pressureFunctionSpace_ = pressureFunctionSpace;
}

//! set the function space object that discretizes the displacements field variable
template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
void QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::
setDisplacementsFunctionSpace(std::shared_ptr<DisplacementsFunctionSpace> displacementsFunctionSpace)
{
  displacementsFunctionSpace_ = displacementsFunctionSpace;

  // also set the functionSpace_ variable which is from the parent class Data
  this->functionSpace_ = displacementsFunctionSpace;
}


template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
void QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::
print()
{
}

template<typename PressureFunctionSpace, typename DisplacementsFunctionSpace>
typename QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::OutputFieldVariables QuasiStaticHyperelasticity<PressureFunctionSpace,DisplacementsFunctionSpace>::
getOutputFieldVariables()
{
  // these field variables will be written to output files
  return std::tuple_cat(
    std::tuple<std::shared_ptr<DisplacementsFieldVariableType>>(std::make_shared<typename DisplacementsFunctionSpace::GeometryFieldType>(this->displacementsFunctionSpace_->geometryField())), // geometry
    std::tuple<std::shared_ptr<DisplacementsFieldVariableType>>(this->displacements_),              // displacements_
    std::tuple<std::shared_ptr<StressFieldVariableType>>(this->pK2Stress_)         // pK2Stress_
  );
}

} // namespace
