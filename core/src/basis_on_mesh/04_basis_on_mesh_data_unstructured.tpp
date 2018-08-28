#include "basis_on_mesh/04_basis_on_mesh_data_unstructured.h"

#include <cmath>
#include <array>
#include <string>
#include <map>
#include <cassert>

#include "easylogging++.h"
#include "utility/string_utility.h"
#include "utility/math_utility.h"

#include "field_variable/unstructured/exfile_representation.h"
#include "field_variable/unstructured/element_to_dof_mapping.h"

namespace BasisOnMesh
{

using namespace StringUtility;

template<int D,typename BasisFunctionType>
BasisOnMeshDataUnstructured<D,BasisFunctionType>::
BasisOnMeshDataUnstructured(std::shared_ptr<Partition::Manager> partitionManager, PyObject *settings, bool noGeometryField) :
  BasisOnMeshPartition<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType>::BasisOnMeshPartition(partitionManager, settings),
  noGeometryField_(noGeometryField)
{
  LOG(TRACE) << "BasisOnMeshDataUnstructured constructor";

}

template<int D,typename BasisFunctionType>
void BasisOnMeshDataUnstructured<D,BasisFunctionType>::
initialize()
{ 
  if (PythonUtility::hasKey(this->specificSettings_, "exelem"))
  {
    std::string filenameExelem = PythonUtility::getOptionString(this->specificSettings_, "exelem", "input.exelem");
    std::string filenameExnode = PythonUtility::getOptionString(this->specificSettings_, "exnode", "input.exnode");

    // read in exelem file
    this->parseExelemFile(filenameExelem);

    // remap names of field variables if specified in config
    // this assigns the geometry field
    this->remapFieldVariables(this->specificSettings_);

    // create the meshPartition by calling BasisOnMeshPartition::initialize() and then create the partitioned vectors for all field variables, also for geometry field
    this->initializeValuesVector();
    
    // read in exnode file
    this->parseExnodeFile(filenameExnode);

    // eliminate scale factors (not yet tested)
    //this->eliminateScaleFactors();
  }
  else if (PythonUtility::hasKey(this->specificSettings_, "nodePositions"))
  {
    // this creates the geometryField and sets the mesh, also creates the meshPartition by calling BasisOnMeshPartition::initialize();
    this->parseFromSettings(this->specificSettings_);
  }
  else
  {
    LOG(FATAL) << "Could not create UnstructuredDeformable node positions. "
      << "Either specify \"exelem\" and \"exnode\" or \"nodePositions\". ";
  }
}

template<int D,typename BasisFunctionType>
dof_no_t BasisOnMeshDataUnstructured<D,BasisFunctionType>::
getDofNo(element_no_t elementNo, int dofIndex) const
{
  if (!this->geometryField_)
    LOG(FATAL) << "Mesh contains no field variable \"geometry\". Use remap to create one!";

  return this->geometryField_->getDofNo(elementNo, dofIndex);
}

template<int D,typename BasisFunctionType>
node_no_t BasisOnMeshDataUnstructured<D,BasisFunctionType>::
getNodeNo(element_no_t elementNo, int nodeIndex) const
{
  return this->elementToNodeMapping_->getElement(elementNo).nodeGlobalNo[nodeIndex];
}

//! return the global node number of element-local node nodeIndex of element with global no elementNoGlobal
template<int D,typename BasisFunctionType>
global_no_t BasisOnMeshDataUnstructured<D,BasisFunctionType>::
getNodeNoGlobal(global_no_t elementNoGlobal, int nodeIndex) const
{
  return this->getNodeNo(elementNoGlobal, nodeIndex);
}

template<int D,typename BasisFunctionType>
void BasisOnMeshDataUnstructured<D,BasisFunctionType>::
getNodeDofs(node_no_t nodeGlobalNo, std::vector<dof_no_t> &dofGlobalNos) const
{
  if (!this->geometryField_)
    LOG(FATAL) << "Mesh contains no field variable \"geometry\". Use remap to create one!";

  std::vector<dof_no_t> &nodeDofs = this->geometryField_->nodeToDofMapping()->getNodeDofs(nodeGlobalNo);

  dofGlobalNos.reserve(dofGlobalNos.size() + nodeDofs.size());

  for(dof_no_t dof : nodeDofs)
  {
    dofGlobalNos.push_back(dof);
  }
}

template<int D,typename BasisFunctionType>
dof_no_t BasisOnMeshDataUnstructured<D,BasisFunctionType>::
getNodeDofNo(node_no_t nodeGlobalNo, int dofIndex) const 
{
  return this->geometryField_->nodeToDofMapping()->getNodeDofNo(nodeGlobalNo, dofIndex);
}

template<int D,typename BasisFunctionType>
std::shared_ptr<FieldVariable::ElementToNodeMapping> BasisOnMeshDataUnstructured<D,BasisFunctionType>::
elementToNodeMapping()
{
  return elementToNodeMapping_;
}

template<int D,typename BasisFunctionType>
element_no_t BasisOnMeshDataUnstructured<D,BasisFunctionType>::
nElementsLocal() const
{
  assert(geometryField_);
  return this->geometryField_->nElements();
}

template<int D,typename BasisFunctionType>
global_no_t BasisOnMeshDataUnstructured<D,BasisFunctionType>::
nElementsGlobal() const
{
  assert(geometryField_);
  return this->geometryField_->nElements();
}

};  // namespace
