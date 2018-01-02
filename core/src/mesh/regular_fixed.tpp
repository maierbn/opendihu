#include "control/types.h"

#include <array>
#include <utility/python_utility.h>
#include "easylogging++.h"

namespace Mesh
{

template<unsigned long D>
RegularFixed<D>::RegularFixed(PyObject *specificSettings) : Structured<D>(specificSettings)
{
  // get settings values nElements_ and physical Extend
  std::array<double, D> defaultValues;
  defaultValues.fill(1.0);
  
  std::array<double, D> physicalExtend;
  // only get physicalExtend if it is not a 1-node mesh with 0 elements
  if (D > 1 || this->nElements_[0] != 0)
    physicalExtend = PythonUtility::getOptionArray<double, D>(specificSettings, "physicalExtend", 1.0, PythonUtility::Positive);
  else
    physicalExtend[0] = 1.0;
 
  LOG(DEBUG) << "  RegularFixed constructor, D="<< D<<", nElements: ";
  for(auto n : this->nElements_)
  {
     LOG(DEBUG) << "    "<<n;
  }
  LOG(DEBUG) << "  physicalExtend: ";
  for(auto n : physicalExtend)
  {
     LOG(DEBUG) << "    "<<n;
  }
  
  // compute mesh widths from physical extent and number of elements in the coordinate directions
  if (D > 1 || this->nElements_[0] != 0)
  {
    auto nElementsIter = this->nElements_.begin();
    auto physicalExtendIter = physicalExtend.begin();
    for (typename std::array<double, D>::iterator meshWidthIter = meshWidth_.begin(); meshWidthIter != meshWidth_.end(); 
        meshWidthIter++, nElementsIter++, physicalExtendIter++)
    {
      *meshWidthIter = *physicalExtendIter / *nElementsIter;
    }
  }
  else
  {
    // 1D 1-node mesh
    meshWidth_[0] = 1.0;
  }
  
  LOG(DEBUG) << "  meshWidth: ";
  for(auto n : meshWidth_)
  {
     LOG(DEBUG) << "    "<<n;
  }
}
  

template<unsigned long D>
RegularFixed<D>::RegularFixed(std::array<element_idx_t, D> nElements, std::array<double, D> physicalExtent) :
  Structured<D>(nElements), meshWidth_(physicalExtent)
{
  // compute mesh widths from physical extent and number of elements in the coordinate directions
  typename std::array<element_idx_t, D>::iterator nElementsIter = this->nElements_.begin();
  for (typename std::array<double, D>::iterator meshWidthIter = meshWidth_.begin(); meshWidthIter != meshWidth_.end(); 
       meshWidthIter++, nElementsIter++)
  {
    *meshWidthIter = *meshWidthIter / *nElementsIter;
  }
}

template<unsigned long D>
double RegularFixed<D>::meshWidth(int dimension) const
{
  return meshWidth_[dimension];
}
 
};