#include "data_management/finite_element_method/diffusion_tensor_base.h"

namespace Data
{


template<typename FunctionSpaceType>
DiffusionTensorBase<FunctionSpaceType>::
DiffusionTensorBase(PythonConfig specificSettings) :
  specificSettings_(specificSettings)
{
}


}  // namespace
