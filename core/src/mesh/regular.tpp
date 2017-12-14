
#include "mesh/regular.h"

#include "utility/python_utility.h"

namespace Mesh
{

template<unsigned long D>
Regular<D>::Regular(std::array<element_idx_t, D> &nElements) :
  MeshD<D>(NULL), nElements_(nElements)
{
}
  
template<unsigned long D>
Regular<D>::Regular(PyObject *specificSettings) : MeshD<D>(specificSettings)
{
  // get settings values nElements_
  if(specificSettings != NULL)
  {
    this->nElements_ = PythonUtility::getOptionArray<int, D>(specificSettings, "nElements", 10, PythonUtility::NonNegative);
  }
} 
  
template<unsigned long D>
element_idx_t Regular<D>::nElements(int dimension)
{
  return nElements_[dimension];
}

template<unsigned long D>
element_idx_t Regular<D>::nNodes(int dimension)
{
  return nElements_[dimension]+1;
}

};    // namespace