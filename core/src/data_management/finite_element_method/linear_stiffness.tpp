#include "data_management/finite_element_method/linear_stiffness.h"

#include <Python.h>  // has to be the first included header

namespace Data
{

template<typename FunctionSpaceType,int nComponents>
void LinearStiffness<FunctionSpaceType,nComponents>::
initialize()
{
  FiniteElementsBase<FunctionSpaceType,nComponents>::initialize();

  // load K and μ from settings
  bulkModulus_ = this->context_.getPythonConfig().getOptionDouble("bulkModulus", 1.0, PythonUtility::ValidityCriterion::Positive);
  shearModulus_ = this->context_.getPythonConfig().getOptionDouble("shearModulus", 1.0, PythonUtility::ValidityCriterion::Positive);

  activeStress_ = nullptr;
}

//! get the value of the 2nd order stiffness tensor, C_abcd
template<typename FunctionSpaceType,int nComponents>
double LinearStiffness<FunctionSpaceType,nComponents>::
linearStiffness(int a, int b, int c, int d) const
{
  // formula: C_abcd = K δ_ab δ_cd + μ(δ_ac δ_bd + δ_ad δ_bc - 2/3 δ_ab δ_cd)

  double delta_ab = (a == b? 1.0 : 0.0);
  double delta_cd = (c == d? 1.0 : 0.0);
  double delta_ac = (a == c? 1.0 : 0.0);
  double delta_bd = (b == d? 1.0 : 0.0);
  double delta_ad = (a == d? 1.0 : 0.0);
  double delta_bc = (b == c? 1.0 : 0.0);

  return bulkModulus_ * delta_ab * delta_cd + shearModulus_ * (delta_ac * delta_bd + delta_ad * delta_bc - 2./3 * delta_ab * delta_cd);
}

//! set values for active stress
template<typename FunctionSpaceType,int nComponents>
void LinearStiffness<FunctionSpaceType,nComponents>::
setActiveStress(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents*nComponents>> activeStress)
{
  activeStress_ = activeStress;
}

  //! get the active stress DxD tensor (row major)
template<typename FunctionSpaceType,int nComponents>
std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents*nComponents>> LinearStiffness<FunctionSpaceType,nComponents>::
activeStress()
{
  return activeStress_;
}

}  // namespace
