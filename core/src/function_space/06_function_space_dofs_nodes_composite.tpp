#include "function_space/06_function_space_dofs_nodes.h"

#include <Python.h>  // has to be the first included header
#include "easylogging++.h"

#include <cmath>
#include <array>
#include <string>

namespace FunctionSpace
{

// constructor
template<int D,typename BasisFunctionType>
FunctionSpaceDofsNodes<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
FunctionSpaceDofsNodes(std::shared_ptr<Partition::Manager> partitionManager,
                       std::vector<std::shared_ptr<FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>> subFunctionSpaces, bool noGeometryField) :
  FunctionSpaceGeometry<Mesh::CompositeOfDimension<D>,BasisFunctionType>(partitionManager, subFunctionSpaces)
{
  this->noGeometryField_ = noGeometryField;
}

//! constructor with PythonConfig, this is needed such that mesh manager compiles but never called
template<int D,typename BasisFunctionType>
FunctionSpaceDofsNodes<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
FunctionSpaceDofsNodes(std::shared_ptr<Partition::Manager> partitionManager,
                        PythonConfig settings, bool noGeometryField):
  FunctionSpaceGeometry<Mesh::CompositeOfDimension<D>,BasisFunctionType>(partitionManager, std::vector<std::shared_ptr<FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>>())
{
  LOG(FATAL) << "Constructor of composite mesh with PythonConfig is not possible.";
}

//! constructor with additionally given node positions, this is needed such that mesh manager compiles but never called
template<int D,typename BasisFunctionType>
FunctionSpaceDofsNodes<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
FunctionSpaceDofsNodes(std::shared_ptr<Partition::Manager> partitionManager, std::vector<double> &nodePositions,
                        PythonConfig settings, bool noGeometryField) :
  FunctionSpaceGeometry<Mesh::CompositeOfDimension<D>,BasisFunctionType>(partitionManager, std::vector<std::shared_ptr<FunctionSpace<Mesh::StructuredDeformableOfDimension<D>,BasisFunctionType>>>())
{
  LOG(FATAL) << "Constructor of composite mesh with PythonConfig and node positions is not possible.";
}

template<int D,typename BasisFunctionType>
void FunctionSpaceDofsNodes<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
initialize()
{
  // create meshPartition and redistribute elements if necessary, this needs information about mesh size
  FunctionSpacePartition<Mesh::CompositeOfDimension<D>,BasisFunctionType>::initialize();

  // if this mesh does not have a geometry field, do nothing further
  if (this->noGeometryField_)
    return;

  // setup geometry field

  // pass a shared "this" pointer to the geometryField
  // retrieve "this" pointer and convert to downwards pointer of most derived class "FunctionSpace"
  std::shared_ptr<FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>> thisFunctionSpace
    = std::static_pointer_cast<FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>>(this->shared_from_this());

  assert(thisFunctionSpace != nullptr);

  typedef FieldVariable::FieldVariable<FunctionSpace<Mesh::CompositeOfDimension<D>,BasisFunctionType>,3> GeometryFieldType;

  // create empty field variable for geometry field
  std::vector<std::string> componentNames{"x", "y", "z"};
  this->geometryField_ = std::make_shared<GeometryFieldType>(thisFunctionSpace, "geometry", componentNames, true);

  LOG(DEBUG) << "create composite geometry field with " << this->nDofsLocalWithoutGhosts() << " local entries.";

  // assign values of geometry field
  std::vector<Vec3> targetGeometryValues(this->nDofsLocalWithoutGhosts());

  // loop over submeshes
  for (int subMeshNo = 0; subMeshNo < this->subFunctionSpaces_.size(); subMeshNo++)
  {
    std::vector<Vec3> geometryValues;
    this->subFunctionSpaces_[subMeshNo]->geometryField().getValuesWithoutGhosts(geometryValues);

    // loop over dofs
    for (node_no_t nodeNoLocal = 0; nodeNoLocal < this->subFunctionSpaces_[subMeshNo]->nNodesLocalWithoutGhosts(); nodeNoLocal++)
    {
      node_no_t compositeNodeNo = this->meshPartition()->getNodeNoLocalFromSubmesh(subMeshNo, nodeNoLocal);

      if (compositeNodeNo != -1)
      {
        for (int nodalDofNo = 0; nodalDofNo < this->subFunctionSpaces_[subMeshNo]->nDofsPerNode(); nodalDofNo++)
        {
          dof_no_t compositeDofNo = compositeNodeNo * this->subFunctionSpaces_[subMeshNo]->nDofsPerNode() + nodalDofNo;
          dof_no_t subFunctionSpaceDofNo = nodeNoLocal * this->subFunctionSpaces_[subMeshNo]->nDofsPerNode() + nodalDofNo;

          LOG(DEBUG) << "  dof " << compositeDofNo << " <- [" << subFunctionSpaceDofNo << "], mesh " << subMeshNo;

          targetGeometryValues[compositeDofNo] = geometryValues[subFunctionSpaceDofNo];
        }
      }
    }
  }

  this->geometryField_->setValuesWithoutGhosts(targetGeometryValues);
  this->geometryField_->finishGhostManipulation();

  // set initalized_ to true which indicates that initialize has been called
  this->initialized_ = true;
}

template<int D,typename BasisFunctionType>
node_no_t FunctionSpaceDofsNodes<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
nNodesLocalWithGhosts() const
{
  return this->meshPartition()->nNodesLocalWithGhosts();
}

template<int D,typename BasisFunctionType>
node_no_t FunctionSpaceDofsNodes<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
nNodesLocalWithoutGhosts() const
{
  return this->meshPartition()->nNodesLocalWithoutGhosts();
}

template<int D,typename BasisFunctionType>
dof_no_t FunctionSpaceDofsNodes<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
nDofsLocalWithGhosts() const
{
  return this->meshPartition()->nDofsLocalWithGhosts();
}

template<int D,typename BasisFunctionType>
dof_no_t FunctionSpaceDofsNodes<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
nDofsLocalWithoutGhosts() const
{
  return this->meshPartition()->nDofsLocalWithoutGhosts();
}

template<int D,typename BasisFunctionType>
global_no_t FunctionSpaceDofsNodes<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
nNodesGlobal() const
{
  return this->meshPartition()->nNodesGlobal();
}

template<int D,typename BasisFunctionType>
global_no_t FunctionSpaceDofsNodes<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
nDofsGlobal() const
{
  return this->meshPartition()->nDofsGlobal();
}

template<int D,typename BasisFunctionType>
void FunctionSpaceDofsNodes<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
getNodePositions(std::vector<double> &nodes) const
{
  // assert that geometry field variable is set
  assert (this->geometryField_);

  nodes.resize(this->nNodesLocalWithGhosts()*3);

  for (node_no_t nodeGlobalNo = 0; nodeGlobalNo < this->nNodesLocalWithGhosts(); nodeGlobalNo++)
  {
    node_no_t nodeFirstDofGlobalNo = this->geometryField_->nodeToDofMapping()->getNodeDofs(nodeGlobalNo)[0];
    Vec3 position = this->geometryField_->template getValue(nodeFirstDofGlobalNo);
    int index = nodeGlobalNo*3;
    nodes[index+0] = position[0];
    nodes[index+1] = position[1];
    nodes[index+2] = position[2];
  }
}


} // namespace
