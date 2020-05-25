#pragma once

#include <Python.h>  // has to be the first included header

#include "utility/math_utility.h"
#include "control/python_config/python_config.h"

namespace Data
{

template<typename FunctionSpaceType>
class DiffusionTensorBase
{
public:

  //! constructor
  DiffusionTensorBase(PythonConfig specificSettings);

protected:

  PythonConfig specificSettings_;                         //< the python settings from which the diffusion tensors are parsed
  std::shared_ptr<FunctionSpaceType> dataFunctionSpace_;  //< the function space of the data class, this cannot be named functionSpace_ because the Data class has a variable with this name
};

}  // namespace

#include "data_management/finite_element_method/diffusion_tensor_base.tpp"
