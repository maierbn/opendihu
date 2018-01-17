#include "basis_on_mesh/04_basis_on_mesh_nodes.h"

#include "easylogging++.h"
#include <array>

#include "utility/python_utility.h"
#include "control/types.h"
#include "field_variable/field_variable.h"
#include "field_variable/field_variable_regular_fixed.h"

namespace BasisOnMesh
{

template<int D,typename BasisFunctionType>
BasisOnMeshNodes<Mesh::RegularFixed<D>,BasisFunctionType>::
BasisOnMeshNodes(PyObject *specificSettings) :
  BasisOnMeshDofs<Mesh::RegularFixed<D>,BasisFunctionType>::BasisOnMeshDofs(specificSettings)
{
  std::array<double,D> meshWidths;
 
  // get settings values nElements_ and physical Extend
  std::array<double, D> defaultValues;
  defaultValues.fill(1.0);
  
  std::array<double, D> physicalExtend;
  // only get physicalExtend if it is not a 1-node mesh with 0 elements
  if (D > 1 || this->nElements_[0] != 0)
    physicalExtend = PythonUtility::getOptionArray<double, D>(specificSettings, "physicalExtend", 1.0, PythonUtility::Positive);
  else
    physicalExtend[0] = 1.0;
 
  
  // compute mesh widths from physical extent and number of elements in the coordinate directions
  if (D > 1 || this->nElements_[0] != 0)
  {
    auto nElementsIter = this->nElements_.begin();
    auto physicalExtendIter = physicalExtend.begin();
    for (typename std::array<double,D>::iterator meshWidthIter = meshWidths.begin(); meshWidthIter != meshWidths.end(); 
        meshWidthIter++, nElementsIter++, physicalExtendIter++)
    {
      *meshWidthIter = *physicalExtendIter / *nElementsIter;
    }
  }
  else
  {
    // 1D 1-node mesh
    meshWidths[0] = 1.0;
  }
  
  LOG(DEBUG) << "  BasisOnMeshNodes Mesh::RegularFixed constructor, D="<< D<<", nElements: "<<this->nElements_;
  LOG(DEBUG) << "  physicalExtend: " << physicalExtend;
  LOG(DEBUG) << "  meshWidth: " << meshWidths;
  
  LOG(DEBUG) << "   create geometry field ";
  this->geometry_ = std::make_unique<FieldVariableType>(); 
  
  // setup geometry field
  this->geometry_->setMeshWidth(meshWidths);
  std::vector<std::string> componentNames{"x","y","z"};
  bool isGeometryField = true;
  Vec values;
  node_idx_t nDofs = this->nDofs();
  int nEntries = nDofs * 3;   // 3 components (x,y,z) for every dof
  this->geometry_->set("geometry", componentNames, this->nElements_, nEntries, isGeometryField, values);
}
 
template<int D,typename BasisFunctionType>
BasisOnMeshNodes<Mesh::RegularFixed<D>,BasisFunctionType>::
BasisOnMeshNodes(std::array<element_idx_t, D> nElements, std::array<double, D> physicalExtent) :
  BasisOnMeshDofs<Mesh::RegularFixed<D>,BasisFunctionType>::BasisOnMeshDofs(nullptr), meshWidth_(physicalExtent)
{
  // compute mesh widths from physical extent and number of elements in the coordinate directions
  typename std::array<element_idx_t, D>::iterator nElementsIter = this->nElements_.begin();
  for (typename std::array<double, D>::iterator meshWidthIter = this->meshWidth_.begin(); meshWidthIter != this->meshWidth_.end(); 
       meshWidthIter++, nElementsIter++)
  {
    *meshWidthIter = *meshWidthIter / *nElementsIter;
  }
}

template<int D,typename BasisFunctionType>
void BasisOnMeshNodes<Mesh::RegularFixed<D>,BasisFunctionType>::
initialize()
{
  LOG(DEBUG) << "   retrieve this pointer ";
  std::shared_ptr<BasisOnMeshNodes<Mesh::RegularFixed<D>,BasisFunctionType>> ptr = this->shared_from_this();
  
  assert(ptr != nullptr);
  
  LOG(DEBUG) << "   cast this pointer ";
  std::shared_ptr<BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>> self = std::static_pointer_cast<BasisOnMesh<Mesh::RegularFixed<D>,BasisFunctionType>>(ptr);
  
  assert(self != nullptr);
  this->geometry_->setMesh(self);
}

template<int D,typename BasisFunctionType>
node_idx_t BasisOnMeshNodes<Mesh::RegularFixed<D>,BasisFunctionType>::
nNodes() const
{
  int result = 1;
  for (int i=0; i<D; i++)
    result *= nNodes(i);
  return result;
}

template<int D,typename BasisFunctionType>
node_idx_t BasisOnMeshNodes<Mesh::RegularFixed<D>,BasisFunctionType>::
nDofs() const
{
  return nNodes() * this->nDofsPerNode();
}

template<int D,typename BasisFunctionType>
node_idx_t BasisOnMeshNodes<Mesh::RegularFixed<D>,BasisFunctionType>::
nNodes(int dimension) const
{
  return this->nElements(dimension) * BasisOnMeshBaseDim<1,BasisFunctionType>::averageNNodesPerElement() + 1;
}

template<int D,typename BasisFunctionType>
double BasisOnMeshNodes<Mesh::RegularFixed<D>,BasisFunctionType>::
meshWidth(int dimension) const
{
  return this->geometry_->meshWidth(dimension);
}
 
template<int D,typename BasisFunctionType>
Vec3 BasisOnMeshNodes<Mesh::RegularFixed<D>,BasisFunctionType>::
getGeometry(node_idx_t dofGlobalNo) const
{
  Vec3 result = geometry_.getValue<3>(dofGlobalNo);
  return result;
}  
  
template<int D,typename BasisFunctionType>
void BasisOnMeshNodes<Mesh::RegularFixed<D>,BasisFunctionType>::
getElementGeometry(element_idx_t elementNo, std::array<Vec3, BasisOnMeshBaseDim<D,BasisFunctionType>::nDofsPerElement()> &values)
{
  const int nDofsPerElement = BasisOnMeshBaseDim<D,BasisFunctionType>::nDofsPerElement();
  geometry_.getElementValues<nDofsPerElement,3>(elementNo, values);
}

template<int D,typename BasisFunctionType>
typename BasisOnMeshNodes<Mesh::RegularFixed<D>,BasisFunctionType>::FieldVariableType &BasisOnMeshNodes<Mesh::RegularFixed<D>,BasisFunctionType>::
geometryField()
{
  if (this->geometry_ == nullptr)
    LOG(ERROR) << "Geometry field is not yet set.";
  return *this->geometry_;
}

template<int D,typename BasisFunctionType>
void BasisOnMeshNodes<Mesh::RegularFixed<D>,BasisFunctionType>::
getNodePositions(std::vector<double> &nodes) const
{
  nodes.resize(this->nNodes()*3);
 
  for (int nodeGlobalNo = 0; nodeGlobalNo < this->nNodes(); nodeGlobalNo++)
  {
    int firstNodeDofGlobalNo = nodeGlobalNo*this->nDofsPerNode();
    
    int index = nodeGlobalNo*3;
    Vec3 position = this->geometry_->template getValue<3>(firstNodeDofGlobalNo);
    nodes[index+0] = position[0];
    nodes[index+1] = position[1];
    nodes[index+2] = position[2];
  }
}


};  // namespace