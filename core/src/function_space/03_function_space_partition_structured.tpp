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
    return;
  
  // Creation of the partitioning is only possible after the number of elements is known.
  // Because this may need file I/O (e.g. reading from exfiles)
 
  VLOG(1) << "FunctionSpacePartition<Structured>::initialize()";
  
  // create partitioning
  assert(this->partitionManager_ != nullptr);
  
  // get if the mesh information in config specifies local or global domain
  bool inputMeshIsGlobal = PythonUtility::getOptionBool(this->specificSettings_, "inputMeshIsGlobal", true);
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
