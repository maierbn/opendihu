#include "data_management/finite_elements.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <memory>

#include "easylogging++.h"

#include "utility/python_utility.h"
#include "control/dihu_context.h"
#include "utility/petsc_utility.h"
#include "function_space/function_space.h"
#include "mesh/unstructured_deformable.h"
#include "basis_function/hermite.h"
#include "partition/partitioned_petsc_mat.h"

namespace Data
{

//! constructor
template<typename FunctionSpaceType>
FiniteElements<FunctionSpaceType,Equation::Dynamic::AnisotropicDiffusion>::
FiniteElements(DihuContext context) :
  FiniteElementsBase<FunctionSpaceType>(context),
  DiffusionTensorConstant<FunctionSpaceType>(context.getPythonConfig())
{
}

//! constructor
template<typename FunctionSpaceType>
FiniteElements<FunctionSpaceType,Equation::Dynamic::DirectionalDiffusion>::
FiniteElements(DihuContext context) :
  FiniteElementsBase<FunctionSpaceType>(context),
  DiffusionTensorFieldVariable<FunctionSpaceType>(context.getPythonConfig())
{
}

template<typename FunctionSpaceType>
void FiniteElements<FunctionSpaceType,Equation::Dynamic::AnisotropicDiffusion>::
initialize()
{
  FiniteElementsBase<FunctionSpaceType>::initialize();

  // set up diffusion tensor if there is any
  DiffusionTensorConstant<FunctionSpaceType>::initialize();
}

template<typename FunctionSpaceType>
void FiniteElements<FunctionSpaceType,Equation::Dynamic::DirectionalDiffusion>::
initialize(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> direction, bool useAdditionalDiffusionTensor)
{
  FiniteElementsBase<FunctionSpaceType>::initialize();

  // set up diffusion tensor, initialize with given direction field
  DiffusionTensorFieldVariable<FunctionSpaceType>::initialize(direction, useAdditionalDiffusionTensor);
}


} // namespace Data
