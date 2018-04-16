#include "basis_on_mesh/05_basis_on_mesh_dofs_nodes.h"

#include <Python.h>  // has to be the first included header
#include "easylogging++.h"
#include <array>

#include "utility/python_utility.h"
#include "control/types.h"
#include "utility/vector_operators.h"
#include "field_variable/field_variable.h"
#include "field_variable/00_field_variable_base.h"

namespace BasisOnMesh
{

template<int D,typename BasisFunctionType>
BasisOnMeshDofsNodes<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>::
BasisOnMeshDofsNodes(PyObject *specificSettings) :
  BasisOnMeshGeometry<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>::BasisOnMeshGeometry(specificSettings)
{
  this->meshWidth_ = 0;
 
  // get settings values nElements_ and physical Extend
  std::array<double, D> defaultValues;
  defaultValues.fill(1.0);
  
  std::array<double, D> physicalExtent;
  // only get physicalExtent if it is not a 1-node mesh with 0 elements
  if (D > 1 || this->nElementsPerCoordinateDirection_[0] != 0)
    physicalExtent = PythonUtility::getOptionArray<double, D>(specificSettings, "physicalExtent", 1.0, PythonUtility::Positive);
  else
    physicalExtent[0] = 1.0;
 
  
  // compute mesh width from physical extent and number of elements in the coordinate directions
  if (D > 1 || this->nElementsPerCoordinateDirection_[0] != 0)
  {
    auto nElementsIter = this->nElementsPerCoordinateDirection_.begin();
    auto physicalExtentIter = physicalExtent.begin();
    for (; physicalExtentIter != physicalExtent.end(); nElementsIter++, physicalExtentIter++)
    {
      double meshWidthCurrentDirection = *physicalExtentIter / *nElementsIter;
      if (this->meshWidth_ == 0)
      {
        this->meshWidth_ = meshWidthCurrentDirection;
      }
      else if(fabs(this->meshWidth_ - meshWidthCurrentDirection) > 1e-14)
      {
        LOG(ERROR) << "Mesh has no uniform mesh width, use a structured deformable grid instead.";
        LOG(ERROR) << "mesh width: " << this->meshWidth_ << ", other: " << *physicalExtentIter << "/" << *nElementsIter << "=" << meshWidthCurrentDirection;
      }
    }
  }
  else
  {
    // 1D 1-node mesh
    this->meshWidth_ = 1.0;
  }
  
  LOG(DEBUG) << "  BasisOnMeshDofsNodes Mesh::RegularFixed constructor, D="<< D<<", nElements: "<<this->nElementsPerCoordinateDirection_;
  LOG(DEBUG) << "  physicalExtent: " << physicalExtent;
  LOG(DEBUG) << "  meshWidth: " << this->meshWidth_;
  
  LOG(DEBUG) << "   create geometry field ";
  
  setupGeometryField();
}

template<int D,typename BasisFunctionType>
BasisOnMeshDofsNodes<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>::
BasisOnMeshDofsNodes(std::array<element_no_t, D> nElements, std::array<double, D> physicalExtent) :
  BasisOnMeshGeometry<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>::BasisOnMeshGeometry(nullptr)
{
  // compute mesh width from physical extent and number of elements in the coordinate directions
  this->meshWidth_ = 0;
  typename std::array<element_no_t, D>::iterator nElementsIter = this->nElementsPerCoordinateDirection_.begin();
  for (typename std::array<double, D>::iterator physicalExtentIter = physicalExtent.begin(); physicalExtentIter != physicalExtent.end(); 
       physicalExtentIter++, nElementsIter++)
  {
    double meshWidthCurrentDirection = *physicalExtentIter / *nElementsIter;
    if (this->meshWidth_ == 0)
    {
      this->meshWidth_ = meshWidthCurrentDirection;
    }
    else if(fabs(this->meshWidth_ - meshWidthCurrentDirection) > 1e-14)
    {
      LOG(ERROR) << "Mesh has no uniform mesh width, use a structured deformable grid instead.";
      LOG(ERROR) << "mesh width: " << this->meshWidth_ << ", other: " << *physicalExtentIter << "/" << *nElementsIter << "=" << meshWidthCurrentDirection;
    }
  }
  
  setupGeometryField();
}

template<int D,typename BasisFunctionType>
void BasisOnMeshDofsNodes<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>::
setupGeometryField()
{
  this->geometryField_ = std::make_unique<GeometryFieldType>();
  
  // setup geometry field
  this->geometryField_->setMeshWidth(this->meshWidth_);
  std::vector<std::string> componentNames{"x","y","z"};
  bool isGeometryField = true;
  Vec values;
  dof_no_t nDofs = this->nDofs();
  std::size_t nEntries = nDofs * 3;   // 3 components (x,y,z) for every dof
  this->geometryField_->set("geometry", componentNames, this->nElementsPerCoordinateDirection_, nEntries, isGeometryField, values);
}


template<int D,typename BasisFunctionType>
node_no_t BasisOnMeshDofsNodes<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>::
nNodes() const
{
  int result = 1;
  for (int i=0; i<D; i++)
    result *= nNodes(i);
  return result;
}

template<int D,typename BasisFunctionType>
dof_no_t BasisOnMeshDofsNodes<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>::
nDofs() const
{
  return nNodes() * this->nDofsPerNode();
}

template<int D,typename BasisFunctionType>
node_no_t BasisOnMeshDofsNodes<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>::
nNodes(int dimension) const
{
  //LOG(DEBUG) << "nNodes (" << dimension << "): " << this->nElementsPerCoordinateDirection(dimension) << "*" << BasisOnMeshBaseDim<1,BasisFunctionType>::averageNNodesPerElement() << "+1";
  return this->nElementsPerCoordinateDirection(dimension) * BasisOnMeshBaseDim<1,BasisFunctionType>::averageNNodesPerElement() + 1;
}

template<int D,typename BasisFunctionType>
double BasisOnMeshDofsNodes<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>::
meshWidth() const
{
  return this->geometryField_->meshWidth();
}

template<int D,typename BasisFunctionType>
void BasisOnMeshDofsNodes<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>::
getNodePositions(std::vector<double> &nodes) const
{
  nodes.resize(this->nNodes()*3);
 
  for (node_no_t nodeGlobalNo = 0; nodeGlobalNo < this->nNodes(); nodeGlobalNo++)
  {
    dof_no_t firstNodeDofGlobalNo = nodeGlobalNo*this->nDofsPerNode();
    
    std::size_t index = nodeGlobalNo*3;
    Vec3 position = this->geometryField_->getValue(firstNodeDofGlobalNo);
    nodes[index+0] = position[0];
    nodes[index+1] = position[1];
    nodes[index+2] = position[2];
  }
}


};  // namespace