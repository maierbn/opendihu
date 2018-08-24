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
 
  // This initialize() method is called from parseFromSettings(), 04_basis_on_mesh_data_unstructured_exfile_io.tpp
  
  VLOG(1) << "BasisOnMeshPartition<Unstructured>::initialize()";
  
  // create partitioning
  assert(this->partitionManager_ != nullptr);
  
  assert(this->nElementsGlobal() != 0);
  assert(this->nNodesGlobal() != 0);
  assert(this->nDofsGlobal() != 0);
  
  this->meshPartition_ = this->partitionManager_->template createPartitioningUnstructured<BasisOnMesh<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>(
    this->nElementsGlobal(), this->nNodesGlobal(), this->nDofsGlobal());
}

};  // namespace