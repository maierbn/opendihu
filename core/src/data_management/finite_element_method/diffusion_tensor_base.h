#pragma once

#include <Python.h>  // has to be the first included header

#include "utility/math_utility.h"

namespace Data
{

template<typename FunctionSpaceType>
class DiffusionTensorBase
{
public:

  //! constructor
  DiffusionTensorBase(PyObject *specificSettings);

protected:

  // parse a diffusion tensor from python config, at the given config
  MathUtility::Matrix<FunctionSpaceType::dim(),FunctionSpaceType::dim()> parseDiffusionTensor(std::string settingsKey);

  PyObject *specificSettings_;    ///< the python settings from which the diffusion tensors are parsed
};

}  // namespace

#include "data_management/finite_element_method/diffusion_tensor_base.tpp"
