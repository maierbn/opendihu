#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>
#include <tuple>

#include "data_management/data.h"
#include "data_management/finite_elements_base.h"
#include "data_management/diffusion_tensor_field_variable.h"
#include "data_management/diffusion_tensor_constant.h"
#include "equation/diffusion.h"

namespace Data
{

template<typename FunctionSpaceType,typename Term,typename = Term,typename = typename FunctionSpaceType::BasisFunction>
class FiniteElements :
  public FiniteElementsBase<FunctionSpaceType>,
  public DiffusionTensorConstant<FunctionSpaceType::dim()>
{
public:
  //! constructor from base class
  using FiniteElementsBase<FunctionSpaceType>::FiniteElementsBase;

  // !intialize base class and diffusion tensor
  virtual void initialize();
};

/** for directional diffusion use the diffusion tensor that depends upon a direction field
 */
template<typename FunctionSpaceType>
class FiniteElements<FunctionSpaceType,Equation::Dynamic::DirectionalDiffusion> :
  public FiniteElementsBase<FunctionSpaceType>,
  public DiffusionTensorFieldVariable<FunctionSpaceType>
{
public:
  //! constructor from base class
  using FiniteElementsBase<FunctionSpaceType>::FiniteElementsBase;

  //! dummy method
  virtual void initialize(){};

  // !intialize base class and diffusion tensor which needs the direction field and the number of compartments in the multidomain context
  virtual void initialize(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> direction, int multidomainNCompartments=0);
};



/*
#include "equation/type_traits.h"

template<typename FunctionSpaceType,typename Term>
class FiniteElements<
  FunctionSpaceType,
  Term,
  Equation::isSolidMechanics<Term>,
  BasisFunction::isNotMixed<typename FunctionSpaceType::BasisFunction>
> :
  public Data<FunctionSpaceType>
{
public:
  void tangentStiffnessMatrix();
};
*/
}  // namespace

#include "data_management/finite_elements.tpp"
#include "data_management/finite_elements_mixed.h"
#include "data_management/finite_elements_solid_mechanics.h"
