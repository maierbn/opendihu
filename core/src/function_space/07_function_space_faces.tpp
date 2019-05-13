#include "function_space/07_function_space_faces.h"

#include <Python.h>  // has to be the first included header
#include <array>
#include "easylogging++.h"

namespace FunctionSpace
{

//! get all dof indices of a face for 1D
template<typename MeshType,typename BasisFunctionType>
void FunctionSpaceFaces<MeshType,BasisFunctionType,Mesh::isDim<1,MeshType>> ::
getFaceDofs(Mesh::face_t face, std::array<dof_no_t,FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerNode()> &dofIndices)
{
  assert(face <= Mesh::face0Plus);
  const int nDofsPerElement = FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerElement();
  const int nDofsPerNode = FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerNode();

  switch(face)
  {
  case Mesh::face0Minus:
    for (int i = 0; i < nDofsPerNode; i++)
    {
      dofIndices[i] = i;
    }
    break;

  case Mesh::face0Plus:
    for (int i = 0; i < nDofsPerNode; i++)
    {
      dofIndices[i] = nDofsPerElement - nDofsPerNode + i;
    }
    break;

  default:
    assert(false);
    break;
  }
}

//! get all dof indices of a face for 2D
template<typename MeshType,typename BasisFunctionType>
void FunctionSpaceFaces<MeshType,BasisFunctionType,Mesh::isDim<2,MeshType>> ::
getFaceDofs(Mesh::face_t face, std::array<dof_no_t,FunctionSpaceBaseDim<1,BasisFunctionType>::nNodesPerElement()*FunctionSpaceBaseDim<MeshType::dim(),BasisFunctionType>::nDofsPerNode()> &dofIndices)
{
  assert(face <= Mesh::face1Plus);

  //   1+
  //   _
  //0-|_|0+
  //   1-
  //
  dof_no_t currentDofIndex = 0;
  const int nNodesPerElement1D = FunctionSpaceBaseDim<1,BasisFunctionType>::nNodesPerElement();
  const int nDofsPerNode = FunctionSpaceBaseDim<MeshType::dim(),BasisFunctionType>::nDofsPerNode();
  const int nDofsPer1DLine = FunctionSpaceBaseDim<1,BasisFunctionType>::nNodesPerElement() * nDofsPerNode;
  //const int nDofsPerNode1D = FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerNode();
  switch(face)
  {
  case Mesh::face0Minus:

    currentDofIndex = 0;
    for (int i = 0; i < nNodesPerElement1D; i++)
    {
      for (int nodalDofIndex = 0; nodalDofIndex < nDofsPerNode; nodalDofIndex++)
      {
        dofIndices[i*nDofsPerNode + nodalDofIndex] = currentDofIndex + nodalDofIndex;
      }
      currentDofIndex += nDofsPer1DLine;
    }
    break;

  case Mesh::face0Plus:

    currentDofIndex = nDofsPer1DLine - nDofsPerNode;
    for (int i = 0; i < nNodesPerElement1D; i++)
    {
      for (int nodalDofIndex = 0; nodalDofIndex < nDofsPerNode; nodalDofIndex++)
      {
        dofIndices[i*nDofsPerNode + nodalDofIndex] = currentDofIndex + nodalDofIndex;
      }
      currentDofIndex += nDofsPer1DLine;
    }
    break;

  case Mesh::face1Minus:

    currentDofIndex = 0;
    for (int i = 0; i < nNodesPerElement1D; i++)
    {
      for (int nodalDofIndex = 0; nodalDofIndex < nDofsPerNode; nodalDofIndex++)
      {
        dofIndices[i*nDofsPerNode + nodalDofIndex] = currentDofIndex + nodalDofIndex;
      }
      currentDofIndex += nDofsPerNode;
    }
    break;

  case Mesh::face1Plus:

    currentDofIndex = FunctionSpaceBaseDim<2,BasisFunctionType>::nDofsPerElement() - nDofsPer1DLine;
    for (int i = 0; i < nNodesPerElement1D; i++)
    {
      for (int nodalDofIndex = 0; nodalDofIndex < nDofsPerNode; nodalDofIndex++)
      {
        dofIndices[i*nDofsPerNode + nodalDofIndex] = currentDofIndex + nodalDofIndex;
      }
      currentDofIndex += nDofsPerNode;
    }
    break;
  default:
    assert(false);
    break;
  }
}

//! get all dof indices of a face for 3D
template<typename MeshType,typename BasisFunctionType>
void FunctionSpaceFaces<MeshType,BasisFunctionType,Mesh::isDim<3,MeshType>> ::
getFaceDofs(Mesh::face_t face, std::array<dof_no_t,FunctionSpaceBaseDim<2,BasisFunctionType>::nNodesPerElement()*FunctionSpaceBaseDim<MeshType::dim(),BasisFunctionType>::nDofsPerNode()> &dofIndices)
{
  // bottom layer
  //   1+
  //   _
  //0-|_|0+  (2-)
  //   1-
  // ------------
  // top layer
  //   1+
  //   _
  //0-|_|0+  (2+)
  //   1-
  //
  dof_no_t currentDofIndex = 0;
  dof_no_t rowStartIndex = 0;
  int index = 0;

  const int nNodesPerElement2D = FunctionSpaceBaseDim<2,BasisFunctionType>::nNodesPerElement();
  const int nDofsPerNode = FunctionSpaceBaseDim<MeshType::dim(),BasisFunctionType>::nDofsPerNode();
  const int nNodesPer1DLine = FunctionSpaceBaseDim<1,BasisFunctionType>::nNodesPerElement();
  const int nDofsPer1DLine = nNodesPer1DLine * nDofsPerNode;
  const int nDofsPer2DPlane = nNodesPerElement2D * nDofsPerNode;

  LOG(DEBUG) << "getFaceDofs, MeshType: " << StringUtility::demangle(typeid(MeshType).name()) << ", BasisFunctionType: " << StringUtility::demangle(typeid(BasisFunctionType).name())
    << ", nDofsPer2DPlane: " << nDofsPer2DPlane;

  switch(face)
  {
  case Mesh::face0Minus:    // left

    currentDofIndex = 0;
    for (int i = 0; i < nNodesPerElement2D; i++)
    {
      for (int nodalDofIndex = 0; nodalDofIndex < nDofsPerNode; nodalDofIndex++)
      {
        assert(i*nDofsPerNode + nodalDofIndex < nDofsPer2DPlane);
        dofIndices[i*nDofsPerNode + nodalDofIndex] = currentDofIndex + nodalDofIndex;
      }
      currentDofIndex += nDofsPer1DLine;
    }
    break;

  case Mesh::face0Plus:   // right

    currentDofIndex = nDofsPer1DLine - nDofsPerNode;
    for (int i = 0; i < nNodesPerElement2D; i++)
    {
      for (int nodalDofIndex = 0; nodalDofIndex < nDofsPerNode; nodalDofIndex++)
      {
        assert(i*nDofsPerNode + nodalDofIndex < nDofsPer2DPlane);
        dofIndices[i*nDofsPerNode + nodalDofIndex] = currentDofIndex + nodalDofIndex;
      }
      currentDofIndex += nDofsPer1DLine;
    }
    break;

  case Mesh::face1Minus:      // front

    currentDofIndex = 0;
    rowStartIndex = 0;
    index = 0;
    for (int i = 0; i < nNodesPer1DLine; i++)
    {
      currentDofIndex = rowStartIndex;
      for (int j = 0; j < nNodesPer1DLine; j++)
      {
        for (int nodalDofIndex = 0; nodalDofIndex < nDofsPerNode; nodalDofIndex++, index++)
        {
          assert(index < nDofsPer2DPlane);
          dofIndices[index] = currentDofIndex + nodalDofIndex;
        }
        currentDofIndex += nDofsPerNode;
      }

      rowStartIndex += nDofsPer2DPlane;
    }
    break;

  case Mesh::face1Plus:     // back

    currentDofIndex = 0;
    rowStartIndex = nDofsPer2DPlane - nDofsPer1DLine;
    index = 0;
    for (int i = 0; i < nNodesPer1DLine; i++)
    {
      currentDofIndex = rowStartIndex;
      for (int j = 0; j < nNodesPer1DLine; j++)
      {
        for (int nodalDofIndex = 0; nodalDofIndex < nDofsPerNode; nodalDofIndex++, index++)
        {
          assert(index < nDofsPer2DPlane);
          dofIndices[index] = currentDofIndex + nodalDofIndex;
        }
        currentDofIndex += nDofsPerNode;
      }

      rowStartIndex += nDofsPer2DPlane;
    }
    break;

  case Mesh::face2Minus:    // bottom

    currentDofIndex = 0;
    for (int i = 0; i < nNodesPerElement2D; i++)
    {
      for (int nodalDofIndex = 0; nodalDofIndex < nDofsPerNode; nodalDofIndex++)
      {
        assert(i*nDofsPerNode + nodalDofIndex < nDofsPer2DPlane);
        dofIndices[i*nDofsPerNode + nodalDofIndex] = currentDofIndex + nodalDofIndex;
      }
      currentDofIndex += nDofsPerNode;
    }
    break;

  case Mesh::face2Plus:     // top

    currentDofIndex = FunctionSpaceBaseDim<3,BasisFunctionType>::nDofsPerElement() - nDofsPer2DPlane;
    for (int i = 0; i < nNodesPerElement2D; i++)
    {
      for (int nodalDofIndex = 0; nodalDofIndex < nDofsPerNode; nodalDofIndex++)
      {
        assert(i*nDofsPerNode + nodalDofIndex < nDofsPer2DPlane);
        dofIndices[i*nDofsPerNode + nodalDofIndex] = currentDofIndex + nodalDofIndex;
      }
      currentDofIndex += nDofsPerNode;
    }
    break;

  default:
    assert(false);
    break;
  }
}

} // namespace
