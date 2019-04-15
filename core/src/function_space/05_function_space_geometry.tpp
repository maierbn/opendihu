#include "function_space/05_function_space_geometry.h"

#include "easylogging++.h"

namespace FunctionSpace
{

template<typename MeshType,typename BasisFunctionType,typename DummyForTraits>
Vec3 FunctionSpaceGeometry<MeshType,BasisFunctionType,DummyForTraits>::
getGeometry(node_no_t dofGlobalNo) const
{
  // assert that geometry field variable is set
  assert (this->geometryField_);

  return this->geometryField_->getValue(dofGlobalNo);
}

//! return an array containing all geometry entries for an element
template<typename MeshType,typename BasisFunctionType,typename DummyForTraits>
void FunctionSpaceGeometry<MeshType,BasisFunctionType,DummyForTraits>::
getElementGeometry(element_no_t elementNoLocal, std::array<Vec3, FunctionSpaceBaseDim<MeshType::dim(),BasisFunctionType>::nDofsPerElement()> &values)
{
  // assert that geometry field variable is set
  assert (this->geometryField_);
  assert (elementNoLocal >= 0);
  if (elementNoLocal >= this->nElementsLocal())
    LOG(ERROR) << "FunctionSpace::getElementGeometry elementNoLocal: " << elementNoLocal << ", nElementsLocal: " << this->nElementsLocal();
  assert (elementNoLocal < this->nElementsLocal());

  this->geometryField_->getElementValues(elementNoLocal, values);
}

template<typename MeshType,typename BasisFunctionType,typename DummyForTraits>
void FunctionSpaceGeometry<MeshType,BasisFunctionType,DummyForTraits>::
extractSurfaceGeometry(const std::array<Vec3, FunctionSpaceBaseDim<MeshType::dim(),BasisFunctionType>::nDofsPerElement()> &geometryVolume, Mesh::face_t face,
                       std::array<Vec3, FunctionSpaceBaseDim<MeshType::dim()-1,BasisFunctionType>::nNodesPerElement()> &geometrySurface)
{
  const int D = MeshType::dim();

  const int nNodesSurface = FunctionSpaceBaseDim<MeshType::dim()-1,BasisFunctionType>::nNodesPerElement();
  const int nNodes1D = FunctionSpaceBaseDim<1,BasisFunctionType>::nNodesPerElement();
  const int nDofsPerNode = FunctionSpaceBaseDim<MeshType::dim(),BasisFunctionType>::nDofsPerNode();

  std::array<int,nNodesSurface> nodeNosSurface;   // mapping from volume nodes to surface nodes

  // store the mapping, for the surface nodes, which are the corresponding volume nodes
  if (D == 1)
  {
    assert(false);  // this does not make sense for D == 1
    if (face == Mesh::face_t::face0Minus)
    {
      nodeNosSurface[0] = 0;
    }
    else if (face == Mesh::face_t::face0Plus)
    {
      nodeNosSurface[0] = 1;
    }
  }
  else if (D == 2)
  {
    // volume is 2D, surface is 1D, i.e. line
    if (face == Mesh::face_t::face0Minus)  // left
    {
      for (int i = 0; i < nNodesSurface; i++)
      {
        nodeNosSurface[i] = i*nNodes1D;
      }
    }
    else if (face == Mesh::face_t::face0Plus)  // right
    {
      for (int i = 0; i < nNodesSurface; i++)
      {
        nodeNosSurface[i] = i*nNodes1D + (nNodes1D-1);
      }
    }
    else if (face == Mesh::face_t::face1Minus)   // bottom
    {
      for (int i = 0; i < nNodesSurface; i++)
      {
        nodeNosSurface[i] = i;
      }
    }
    else if (face == Mesh::face_t::face1Plus)   // top
    {
      for (int i = 0; i < nNodesSurface; i++)
      {
        nodeNosSurface[i] = (nNodes1D-1)*nNodes1D + i;
      }
    }
    else assert(false);
  }
  else if (D == 3)
  {
    // volume is 3D, surface is 2D
    assert(nNodesSurface == nNodes1D*nNodes1D);
    if (face == Mesh::face_t::face0Minus)  // left
    {
      for (int i = 0; i < nNodesSurface; i++)
      {
        nodeNosSurface[i] = i*nNodes1D;
      }
    }
    else if (face == Mesh::face_t::face0Plus)  // right
    {
      for (int i = 0; i < nNodesSurface; i++)
      {
        nodeNosSurface[i] = i*nNodes1D + (nNodes1D-1);
      }
    }
    else if (face == Mesh::face_t::face1Minus)   // front
    {
      assert(nNodesSurface == nNodes1D*nNodes1D);
      int index = 0;
      for (int k = 0; k < nNodes1D; k++)
      {
        for (int i = 0; i < nNodes1D; i++, index++)
        {
          nodeNosSurface[index] = k*nNodes1D*nNodes1D + i;
        }
      }
    }
    else if (face == Mesh::face_t::face1Plus)   // back
    {
      int index = 0;
      for (int k = 0; k < nNodes1D; k++)
      {
        for (int i = 0; i < nNodes1D; i++, index++)
        {
          nodeNosSurface[index] = k*nNodes1D*nNodes1D + (nNodes1D-1)*nNodes1D + i;
        }
      }
    }
    else if (face == Mesh::face_t::face2Minus)   // bottom
    {
      for (int i = 0; i < nNodesSurface; i++)
      {
        nodeNosSurface[i] = i;
      }
    }
    else if (face == Mesh::face_t::face2Plus)   // top
    {
      for (int i = 0; i < nNodesSurface; i++)
      {
        nodeNosSurface[i] = (nNodes1D-1)*nNodes1D*nNodes1D + i;
      }
    }
  }

  for (int surfaceNodeIndex = 0; surfaceNodeIndex < nNodesSurface; surfaceNodeIndex++)
  {
    int nodalDofIndex = 0;
    int volumeNodeIndex = nodeNosSurface[surfaceNodeIndex];
    geometrySurface[surfaceNodeIndex] = geometryVolume[volumeNodeIndex*nDofsPerNode + nodalDofIndex];
  }
}


template<typename MeshType,typename BasisFunctionType,typename DummyForTraits>
bool FunctionSpaceGeometry<MeshType,BasisFunctionType,DummyForTraits>::
hasGeometryField()
{
  return this->geometryField_ != nullptr;
}

//! create a non-geometry field field variable with no values being set, with given component names
template<typename MeshType,typename BasisFunctionType,typename DummyForTraits>
typename FunctionSpaceGeometry<MeshType,BasisFunctionType,DummyForTraits>::GeometryFieldType &
FunctionSpaceGeometry<MeshType,BasisFunctionType,DummyForTraits>::
geometryField()
{
  // assert that geometry field variable is set
  assert (this->geometryField_);

  return *this->geometryField_;
}

} // namespace
