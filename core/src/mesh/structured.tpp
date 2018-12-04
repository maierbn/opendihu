
#include "mesh/structured.h"

#include "easylogging++.h"
#include "utility/python_utility.h"
#include "control/python_config.h"

namespace Mesh
{

template<int D>
Structured<D>::Structured(PythonConfig specificSettings) :
  MeshOfDimension<D>(specificSettings), 
  nElementsPerCoordinateDirectionLocal_({0}),
  nElementsPerCoordinateDirectionGlobal_({0})
{
  // get if the mesh information in config specifies local or global domain
  bool inputMeshIsGlobal = specificSettings.getOptionBool("inputMeshIsGlobal", true);
  if (inputMeshIsGlobal)
  {
    // get settings values nElements_
    if (D == 1)
    {
      // for D=1 default value is 0, i.e. if nothing is given, create a 1-dof degenerate mesh
      this->nElementsPerCoordinateDirectionGlobal_ = specificSettings.getOptionArray<global_no_t, D>("nElements", 0, PythonUtility::NonNegative);
    }
    else
    {
      // for D>1 default values is 1, i.e. one element layer in the degenerate dimensions that were not specified
      this->nElementsPerCoordinateDirectionGlobal_ = specificSettings.getOptionArray<global_no_t, D>("nElements", 1, PythonUtility::NonNegative);
    }
    LOG(DEBUG) << "set global number of elements from settings: " << this->nElementsPerCoordinateDirectionGlobal_;
  }
  else 
  {
    // get settings values nElements_
    this->nElementsPerCoordinateDirectionLocal_ = specificSettings.getOptionArray<element_no_t, D>("nElements", 1, PythonUtility::NonNegative);
    LOG(DEBUG) << "set local number of elements from settings: " << this->nElementsPerCoordinateDirectionLocal_;
    
    if (!specificSettings.hasKey("nRanks"))
    {
      LOG(ERROR) << "Config does not specifiy \"nRanks\", a list of the number of ranks in each coordinate direction.";
    }
    this->nRanks_ = specificSettings.getOptionArray<int, D>("nRanks", 1, PythonUtility::Positive);
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
