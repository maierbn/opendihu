#pragma once

#include <Python.h>  // has to be the first included header

namespace Data
{

/** class that holds the linear stiffness parameters
 */
template<typename FunctionSpaceType,int nComponents>
class LinearStiffness :
  public FiniteElementsBase<FunctionSpaceType,nComponents>
{
public:
  //! constructor
  using FiniteElementsBase<FunctionSpaceType,nComponents>::FiniteElementsBase;

  using FiniteElementsBase<FunctionSpaceType,nComponents>::FieldVariablesForOutputWriter;

  //! initialize stifness parameters, then call the initialize method of the base class
  virtual void initialize();

  //! get the value of the 2nd order stiffness tensor, C_abcd
  double linearStiffness(int a, int b, int c, int d) const;

  //! set values for active stress
  void setActiveStress(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents*nComponents>> activeStress);

  //! get the active stress DxD tensor (row major)
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents*nComponents>> activeStress();

  //! set the pointer to rightHandSideActive, the actual variable is stored in the quasi_static_linear_elasticity class
  void setRightHandSideActive(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> rightHandSideActive);

  //! get a debugging field variable that contains the rhs contribution from the active stress
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> rightHandSideActive();

protected:
  double bulkModulus_;   //< material parameter bulk modulus, symbol K
  double shearModulus_;  //< material parameter shear modulus, symbol μ

  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents*nComponents>> activeStress_;  //< active stress field variable, this is a DxD tensor, stored row-major
  std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,nComponents>> rightHandSideActive_;  //< debugging field variable that contains the rhs contribution from the active stress, the actual object is created in quasi_static_linear_elasticity and passed to this variable by setRightHandSideActive
};

}  // namespace

#include "data_management/finite_element_method/linear_stiffness.tpp"
