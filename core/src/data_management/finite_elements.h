#pragma once

#include <Python.h>  // has to be the first included header
#include <petscmat.h>
#include <memory>
#include <tuple>

#include "data_management/data.h"
#include "data_management/diffusion_tensor_field_variable.h"
#include "data_management/diffusion_tensor_constant.h"
#include "control/types.h"
#include "mesh/mesh.h"
#include "field_variable/field_variable.h"
#include "function_space/mixed_function_space.h"
#include "partition/partitioned_petsc_vec.h"
#include "partition/partitioned_petsc_mat.h"

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
class FiniteElements<FunctionSpaceType,Term::DirectionalDiffusion> :
  public FiniteElementsBase<FunctionSpaceType>,
  public DiffusionTensorFieldVariable<FunctionSpaceType::dim()>
{
public:
  //! constructor from base class
  using FiniteElementsBase<FunctionSpaceType>::FiniteElementsBase;

  // !intialize base class and diffusion tensor
  virtual void initialize(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> direction);
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
