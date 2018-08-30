
#include "mesh/structured.h"

#include "easylogging++.h"
#include "utility/python_utility.h"

namespace Mesh
{

template<int D>
Structured<D>::Structured(PyObject *specificSettings) : 
  MeshOfDimension<D>(specificSettings), 
  nElementsPerCoordinateDirectionLocal_({0}),
  nElementsPerCoordinateDirectionGlobal_({0})
{
  // get if the mesh information in config specifies local or global domain
  bool inputMeshIsGlobal = PythonUtility::getOptionBool(specificSettings, "inputMeshIsGlobal", true);
  if (inputMeshIsGlobal)
  {
    // get settings values nElements_
    this->nElementsPerCoordinateDirectionGlobal_ = PythonUtility::getOptionArray<global_no_t, D>(specificSettings, "nElements", 1, PythonUtility::Positive);
    LOG(DEBUG) << "set global number of elements from settings: " << this->nElementsPerCoordinateDirectionGlobal_;
  }
  else 
  {
    // get settings values nElements_
    this->nElementsPerCoordinateDirectionLocal_ = PythonUtility::getOptionArray<element_no_t, D>(specificSettings, "nElements", 1, PythonUtility::NonNegative);
    LOG(DEBUG) << "set local number of elements from settings: " << this->nElementsPerCoordinateDirectionLocal_;
    
    this->nRanks_ = PythonUtility::getOptionArray<int, D>(specificSettings, "nRanks", 1, PythonUtility::Positive);
    LOG(DEBUG) << "set number of ranks in the directions from settings: " << this->nRanks_;
  }
}

template<int D>
element_no_t Structured<D>::
nElementsPerCoordinateDirectionLocal(int dimension) const
{
  if (dimension >= D)
  {
    return 1;
  }
  return nElementsPerCoordinateDirectionLocal_[dimension];
}

template<int D>
std::array<element_no_t, D> Structured<D>::
nElementsPerCoordinateDirectionLocal() const
{
  return nElementsPerCoordinateDirectionLocal_;
}

template<int D>
global_no_t Structured<D>::
nElementsPerCoordinateDirectionGlobal(int dimension) const
{
  if (dimension >= D)
  {
    return 1;
  }
  return nElementsPerCoordinateDirectionGlobal_[dimension];
}

template<int D>
std::array<global_no_t, D> Structured<D>::
nElementsPerCoordinateDirectionGlobal() const
{
  return nElementsPerCoordinateDirectionGlobal_;
}

};    // namespace
