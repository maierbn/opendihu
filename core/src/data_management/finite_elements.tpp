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

template<typename FunctionSpaceType,typename Term,typename DummyForTraits,typename DummyForTraits2>
void FiniteElements<FunctionSpaceType,Term,DummyForTraits,DummyForTraits2>::
initialize()
{
  FiniteElementsBase<FunctionSpaceType>::initialize();

  // set up diffusion tensor if there is any
  DiffusionTensorConstant<FunctionSpaceType::dim()>::initialize(this->context_.getPythonConfig());
}

template<typename FunctionSpaceType>
void FiniteElements<FunctionSpaceType,Equation::Dynamic::DirectionalDiffusion>::
initialize(std::shared_ptr<FieldVariable::FieldVariable<FunctionSpaceType,3>> direction, int multidomainNCompartments)
{
  FiniteElementsBase<FunctionSpaceType>::initialize();

  // set up diffusion tensor, initialize with given direction field
  DiffusionTensorFieldVariable<FunctionSpaceType>::initialize(direction, multidomainNCompartments);
}


} // namespace Data
