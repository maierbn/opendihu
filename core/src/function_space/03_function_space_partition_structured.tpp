#include "function_space/03_function_space_partition.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"

namespace FunctionSpace
{
 
template<typename MeshType,typename BasisFunctionType>
void FunctionSpacePartition<MeshType,BasisFunctionType>::
initialize()
{
  // if meshPartition was already created earlier, do nothing
  if (this->meshPartition_) 
  {
    LOG(DEBUG) << "in FunctionSpacePartition<structured>: meshPartition already set";
    return;
  }
  
  // Creation of the partitioning is only possible after the number of elements is known.
  // Because this may need file I/O (e.g. reading from exfiles)
 
  LOG(DEBUG) << "FunctionSpacePartition<Structured>::initialize(), create meshPartition";
  
  // create partitioning
  assert(this->partitionManager_ != nullptr);
  
  // get whether the mesh information in config specifies local or global domain
  bool inputMeshIsGlobal = this->specificSettings_.getOptionBool("inputMeshIsGlobal", true);
  if (inputMeshIsGlobal)
  {
    this->meshPartition_ = this->partitionManager_->template createPartitioningStructuredGlobal<FunctionSpace<MeshType,BasisFunctionType>>(
      this->nElementsPerCoordinateDirectionGlobal_, this->nElementsPerCoordinateDirectionLocal_, this->nRanks_);
  }
  else 
  {
    this->meshPartition_ = this->partitionManager_->template createPartitioningStructuredLocal<FunctionSpace<MeshType,BasisFunctionType>>(
      this->nElementsPerCoordinateDirectionGlobal_, this->nElementsPerCoordinateDirectionLocal_, this->nRanks_);
  }
  assert(this->meshPartition_);
}

};  // namespace
