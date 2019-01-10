#include "function_space/03_function_space_partition.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"

namespace FunctionSpace
{
 
// forward declaration
template<typename MeshType,typename BasisFunctionType>
class FunctionSpace;

 
template<int D,typename BasisFunctionType>
void FunctionSpacePartition<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::
initialize()
{
  // Creation of the partitioning is only possible after the number of elements is known.
  // Because this may need file I/O (e.g. reading from exfiles)
 
  // This initialize() method is called from parseFromSettings(), 04_function_space_data_unstructured_exfile_io.tpp
  
  VLOG(1) << "FunctionSpacePartition<Unstructured>::initialize()";
  
  // create partitioning
  assert(this->partitionManager_ != nullptr);
  
  assert(this->nElementsGlobal() != 0);
  assert(this->nNodesGlobal() != 0);
  assert(this->nDofsGlobal() != 0);
  
  this->meshPartition_ = this->partitionManager_->template createPartitioningUnstructured<FunctionSpace<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>>(
    this->nElementsGlobal(), this->nNodesGlobal(), this->nDofsGlobal());

  // set initalized_ to true which indicates that initialize has been called
  this->initialized_ = true;
}

};  // namespace
