#include "basis_on_mesh/03_basis_on_mesh_partition.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"

namespace BasisOnMesh
{
 
template<typename MeshType,typename BasisFunctionType>
void BasisOnMeshPartition<MeshType,BasisFunctionType>::
initialize()
{
  // Creation of the partitioning is only possible after the number of elements is known.
  // Because this may need file I/O (e.g. reading from exfiles)
 
  // create partitioning
  assert(this->partitionManager_ != nullptr);
  
  // get if the mesh information in config specifies local or global domain
  std::string inputMeshIsGlobal = PythonUtility::getOptionBool(this->specificSettings, "inputMeshIsGlobal", true);
  if (inputMeshIsGlobal)
  {
    meshPartition_ = partitionManager_->createPartitioningStructuredGlobal(this->nElementsPerCoordinateDirectionGlobal_);
  }
  else 
  {
    meshPartition_ = partitionManager_->createPartitioningStructuredLocal(this->nElementsPerCoordinateDirectionLocal_, this->nRanks_);
  }
}

};  // namespace