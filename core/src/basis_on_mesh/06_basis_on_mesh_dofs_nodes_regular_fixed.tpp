#include "basis_on_mesh/06_basis_on_mesh_dofs_nodes.h"

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
BasisOnMeshDofsNodes(std::shared_ptr<Partition::Manager> partitionManager, PyObject *specificSettings) :
  BasisOnMeshDofsNodesStructured<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>::BasisOnMeshDofsNodesStructured(partitionManager, specificSettings)
{
  this->meshWidth_ = 0;

  // get settings values nElements_ and physical Extend
  std::array<double, D> defaultValues;
  defaultValues.fill(1.0);

  std::array<double, D> physicalExtent;
  // only get physicalExtent if it is not a 1-node mesh with 0 elements
  if (D > 1 || this->nElementsPerCoordinateDirectionLocal_[0] != 0)
  {
    if (PythonUtility::hasKey(specificSettings, "physicalExtend"))
      LOG(ERROR) << "You misspelled \"physicalExtent\" as \"physicalExtend\"!";
    physicalExtent = PythonUtility::getOptionArray<double, D>(specificSettings, "physicalExtent", 1.0, PythonUtility::Positive);
  }
  else
    physicalExtent[0] = 1.0;


  // compute mesh width from physical extent and number of elements in the coordinate directions
  // note for quadratic elements the mesh width is the distance between the nodes, not length of elements
  if (D > 1 || this->nElementsPerCoordinateDirectionLocal_[0] != 0)
  {
    auto nElementsIter = this->nElementsPerCoordinateDirectionLocal_.begin();
    auto physicalExtentIter = physicalExtent.begin();
    for (; physicalExtentIter != physicalExtent.end(); nElementsIter++, physicalExtentIter++)
    {
      double meshWidthCurrentDirection = *physicalExtentIter / (*nElementsIter * BasisOnMeshBaseDim<1,BasisFunctionType>::averageNNodesPerElement());
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

  LOG(DEBUG) << "  BasisOnMeshDofsNodes Mesh::RegularFixed constructor, D="<< D<<", nElements: "<<this->nElementsPerCoordinateDirectionLocal_;
  LOG(DEBUG) << "  physicalExtent: " << physicalExtent;
  LOG(DEBUG) << "  meshWidth: " << this->meshWidth_;
}

template<int D,typename BasisFunctionType>
BasisOnMeshDofsNodes<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>::
BasisOnMeshDofsNodes(std::shared_ptr<Partition::Manager> partitionManager, std::array<element_no_t, D> nElements, std::array<double, D> physicalExtent) :
  BasisOnMeshDofsNodesStructured<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>::BasisOnMeshDofsNodesStructured(partitionManager, nullptr)
{
  // compute mesh width from physical extent and number of elements in the coordinate directions
  // note for quadratic elements the mesh width is the distance between the nodes, not length of elements
  this->meshWidth_ = 0;
  typename std::array<element_no_t, D>::iterator nElementsIter = this->nElementsPerCoordinateDirectionLocal_.begin();
  for (typename std::array<double, D>::iterator physicalExtentIter = physicalExtent.begin(); physicalExtentIter != physicalExtent.end();
       physicalExtentIter++, nElementsIter++)
  {
    double meshWidthCurrentDirection = *physicalExtentIter / (*nElementsIter * BasisOnMeshBaseDim<1,BasisFunctionType>::averageNNodesPerElement());
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

template<int D,typename BasisFunctionType>
void BasisOnMeshDofsNodes<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>::
initialize()
{ 
  // create meshPartition and redistribute elements if necessary, this needs information about mesh size
  BasisOnMeshPartition<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>::initialize();
 
  // pass a shared "this" pointer to the geometryField
  if (this->noGeometryField_)
    return;

  // retrieve "this" pointer and convert to downwards pointer of most derived class "BasisOnMesh"
  std::shared_ptr<BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>> thisMesh
    = std::static_pointer_cast<BasisOnMesh<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>>(this->shared_from_this());

  assert(thisMesh != nullptr);
  
  // create empty field variable for geometry field
  std::vector<std::string> componentNames{"x", "y", "z"};
  this->geometryField_ = std::make_unique<GeometryFieldType>(thisMesh, "geometry", componentNames);
  
  // no need to set values of the geometry field, because there is no data explicitly stored
}

template<int D,typename BasisFunctionType>
double BasisOnMeshDofsNodes<Mesh::StructuredRegularFixedOfDimension<D>,BasisFunctionType>::
meshWidth() const
{
  return this->meshWidth_;
}


};  // namespace