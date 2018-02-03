
#include "mesh/structured.h"

#include "easylogging++.h"
#include "utility/python_utility.h"

namespace Mesh
{

template<int D>
Structured<D>::Structured(std::array<element_no_t, D> &nElements) :
  MeshD<D>(NULL), nElements_(nElements)
{
}
  
template<int D>
Structured<D>::Structured(PyObject *specificSettings) : MeshD<D>(specificSettings)
{
  // get settings values nElements_
  this->nElements_ = PythonUtility::getOptionArray<element_no_t, D>(specificSettings, "nElements", 10, PythonUtility::NonNegative);
} 
  
template<int D>
element_no_t Structured<D>::nElements(int dimension) const
{
  if (dimension >= D)
    return 1;
  return nElements_[dimension];
}


};    // namespace