#include "control/types.h"

#include <array>
#include <control/python_utility.h>
#include "easylogging++.h"

namespace Mesh
{


template<unsigned long D>
element_idx_t RegularFixed<D>::nElements(int dimension)
{
  return nElements_[dimension];
}

template<unsigned long D>
RegularFixed<D>::RegularFixed(PyObject *specificSettings) : MeshD<D>(specificSettings)
{
  // get settings values
  std::array<double, D> defaultValues;
  defaultValues.fill(1.0);
  
  std::array<double, D> physicalExtend = PythonUtility::getOptionList<double, D>(specificSettings, "physicalExtend", 
                                                                                 1.0, PythonUtility::Positive);
  nElements_ = PythonUtility::getOptionList<int, D>(specificSettings, "nElements", 10, PythonUtility::Positive);
 
  LOG(DEBUG) << "RegularFixed constructor, D="<< D<<", nElements: ";
  for(auto n : nElements_)
  {
     LOG(DEBUG) << "  "<<n;
  }
  LOG(DEBUG) << "physicalExtend: ";
  for(auto n : physicalExtend)
  {
     LOG(DEBUG) << "  "<<n;
  }
  
  // compute mesh widths from physical extent and number of elements in the coordinate directions
  auto nElementsIter = nElements_.begin();
  auto physicalExtendIter = physicalExtend.begin();
  for (typename std::array<double, D>::iterator meshWidthIter = meshWidth_.begin(); meshWidthIter != meshWidth_.end(); 
       meshWidthIter++, nElementsIter++, physicalExtendIter++)
  {
    *meshWidthIter = *physicalExtendIter / *nElementsIter;
  }
}
  

template<unsigned long D>
RegularFixed<D>::RegularFixed(std::array<element_idx_t, D> nElements, std::array<double, D> physicalExtent) :
  nElements_(nElements), meshWidth_(physicalExtent)
{
  // compute mesh widths from physical extent and number of elements in the coordinate directions
  auto nElementsIter = nElements_.begin();
  for (auto &meshWidthIter = meshWidth_.begin(); meshWidthIter != meshWidth_.end(); meshWidthIter++, nElementsIter++)
  {
    *meshWidthIter = *meshWidthIter / *nElementsIter;
  }
}

template<unsigned long D>
double RegularFixed<D>::meshWidth(int dimension)
{
  return meshWidth_[dimension];
}
 
};