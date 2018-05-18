
#include "mesh/structured.h"

#include "easylogging++.h"
#include "utility/python_utility.h"

namespace Mesh
{

template<int D>
Structured<D>::Structured(std::array<element_no_t, D> &nElementsPerCoordinateDirection) :
  MeshOfDimension<D>(NULL), nElementsPerCoordinateDirection_(nElementsPerCoordinateDirection)
{
}

template<int D>
Structured<D>::Structured(PyObject *specificSettings) : MeshOfDimension<D>(specificSettings)
{
  // get settings values nElements_
  this->nElementsPerCoordinateDirection_ = PythonUtility::getOptionArray<element_no_t, D>(specificSettings, "nElements", 1, PythonUtility::NonNegative);
  LOG(DEBUG) << "set number of elements from settings: " << this->nElementsPerCoordinateDirection_;
}

template<int D>
element_no_t Structured<D>::nElementsPerCoordinateDirection(int dimension) const
{
  if (dimension >= D)
  {
    return 1;
  }
  return nElementsPerCoordinateDirection_[dimension];
}

template<int D>
std::array<element_no_t, D> Structured<D>::
nElementsPerCoordinateDirection() const
{
  return nElementsPerCoordinateDirection_;
}

};    // namespace