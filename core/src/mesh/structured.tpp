
#include "mesh/structured.h"

#include "easylogging++.h"
#include "utility/python_utility.h"

namespace Mesh
{

template<int D>
Structured<D>::Structured(std::array<element_no_t, D> &nElementsPerDimension) :
  MeshD<D>(NULL), nElementsPerDimension_(nElementsPerDimension)
{
}
  
template<int D>
Structured<D>::Structured(PyObject *specificSettings) : MeshD<D>(specificSettings)
{
  // get settings values nElements_
  this->nElementsPerDimension_ = PythonUtility::getOptionArray<element_no_t, D>(specificSettings, "nElements", 1, PythonUtility::NonNegative);
  LOG(DEBUG) << "set number of elements from settings: " << this->nElementsPerDimension_;
} 
  
template<int D>
element_no_t Structured<D>::nElementsPerDimension(int dimension) const
{
  if (dimension >= D)
    return 1;
  return nElementsPerDimension_[dimension];
}

template<int D>
std::array<element_no_t, D> Structured<D>::
nElementsPerDimension() const
{
  return nElementsPerDimension_;
}

};    // namespace