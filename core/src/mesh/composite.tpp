#include "mesh/composite.h"

namespace Mesh
{

template<int D>
CompositeOfDimension<D>::
CompositeOfDimension(PythonConfig specificSettings) :
  MeshOfDimension<D>(specificSettings)
{
}

//! get the total number of elements in the local domain
template<int D>
element_no_t CompositeOfDimension<D>::
nElementsLocal() const
{
  return nElementsLocal_;
}

//! get the total number of elements in the global domain
template<int D>
global_no_t CompositeOfDimension<D>::
nElementsGlobal() const
{
  return nElementsGlobal_;
}

}  // namespace
