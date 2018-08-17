#include "basis_on_mesh/03_basis_on_mesh_partition.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"

namespace BasisOnMesh
{
 
// forward declaration
template<typename MeshType,typename BasisFunctionType>
class BasisOnMesh;

 
template<int D,typename BasisFunctionType>
void BasisOnMeshPartition<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
initialize()
{
  // Creation of the partitioning is only possible after the number of elements is known.
  // Because this may need file I/O (e.g. reading from exfiles)
 
  VLOG(1) << "BasisOnMeshPartition<Unstructured>::initialize()";
  
  // create partitioning
  assert(this->partitionManager_ != nullptr);
  this->meshPartition_ = this->partitionManager_->template createPartitioning<BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>(
    this->nElementsGlobal());
}

};  // namespace