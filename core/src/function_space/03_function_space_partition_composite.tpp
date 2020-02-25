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
void FunctionSpacePartition<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
initialize()
{
  PyObject *pyObject = this->specificSettings_.pyObject();

  int nSubmeshes = 1;
  if (PyList_Check(pyObject))
  {
    nSubmeshes = PyList_Size(pyObject);
  }
  LOG(DEBUG) << "initialize " << nSubmeshes << " submeshes";
  subFunctionSpaces_.resize(nSubmeshes);

  typedef FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType> SubFunctionSpaceType;

  // created sub meshes
  for (int i = 0; i < nSubmeshes; i++)
  {
    subFunctionSpaces_[i] = std::make_shared<SubFunctionSpaceType>(this->partitionManager_, PythonConfig(this->specificSettings_,i));
  }

  // Creation of the partitioning is only possible after the number of elements is known.
  // Because this may need file I/O (e.g. reading from exfiles)
 
  // This initialize() method is called from parseFromSettings(), 04_function_space_data_unstructured_exfile_io.tpp
  
  VLOG(1) << "FunctionSpacePartition<Composite>::initialize()";

  // initialize number of local and global elements
  this->nElementsLocal_ = 0;
  this->nElementsGlobal_ = 0;

  // iterate over submeshes
  for(std::shared_ptr<SubFunctionSpaceType> &subFunctionSpace : subFunctionSpaces_)
  {
    subFunctionSpace->initialize();

    // count number of elements
    this->nElementsLocal_ += subFunctionSpace->nElementsLocal();
    this->nElementsGlobal_ += subFunctionSpace->nElementsGlobal();
  }

  // create partitioning
  assert(this->partitionManager_ != nullptr);
  this->meshPartition_ = this->partitionManager_->template createPartitioningComposite<BasisFunctionType,D>(subFunctionSpaces_);

  assert(this->meshPartition_);

  // set initalized_ to true which indicates that initialize has been called
  this->initialized_ = true;
}

template<int D,typename BasisFunctionType>
const std::vector<std::shared_ptr<FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>> &
FunctionSpacePartition<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
subFunctionSpaces()
{
  return subFunctionSpaces_;
}

} // namespace
