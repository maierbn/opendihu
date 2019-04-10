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
  const int nDofsPer1DLine = FunctionSpaceBaseDim<1,BasisFunctionType>::nNodesPerElement() * nDofsPerNode;
  const int nDofsPer2DPlane = FunctionSpaceBaseDim<2,BasisFunctionType>::nNodesPerElement() * nDofsPerNode;

  switch(face)
  {
  case Mesh::face0Minus:

    currentDofIndex = 0;
    for (int i = 0; i < nNodesPerElement2D; i++)
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
    for (int i = 0; i < nNodesPerElement2D; i++)
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
    rowStartIndex = 0;
    index = 0;
    for (int i = 0; i < nNodesPerElement2D; i++)
    {
      currentDofIndex = rowStartIndex;
      for (int j = 0; j < nDofsPer1DLine; j++)
      {
        for (int nodalDofIndex = 0; nodalDofIndex < nDofsPerNode; nodalDofIndex++)
        {
          dofIndices[index++] = currentDofIndex + nodalDofIndex;
        }
        currentDofIndex += nDofsPerNode;
      }

      rowStartIndex += nDofsPer2DPlane;
    }
    break;

  case Mesh::face1Plus:

    currentDofIndex = 0;
    rowStartIndex = nDofsPer2DPlane - nDofsPer1DLine;
    index = 0;
    for (int i = 0; i < nNodesPerElement2D; i++)
    {
      currentDofIndex = rowStartIndex;
      for (int j = 0; j < nDofsPer1DLine; j++)
      {
        for (int nodalDofIndex = 0; nodalDofIndex < nDofsPerNode; nodalDofIndex++)
        {
          dofIndices[index++] = currentDofIndex + nodalDofIndex;
        }
        currentDofIndex += nDofsPerNode;
      }

      rowStartIndex += nDofsPer2DPlane;
    }
    break;

  case Mesh::face2Minus:

    currentDofIndex = 0;
    for (int i = 0; i < nNodesPerElement2D; i++)
    {
      for (int nodalDofIndex = 0; nodalDofIndex < nDofsPerNode; nodalDofIndex++)
      {
        dofIndices[i*nDofsPerNode + nodalDofIndex] = currentDofIndex + nodalDofIndex;
      }
      currentDofIndex += nDofsPerNode;
    }
    break;

  case Mesh::face2Plus:

    currentDofIndex = FunctionSpaceBaseDim<3,BasisFunctionType>::nDofsPerElement() - nDofsPer2DPlane;
    for (int i = 0; i < nNodesPerElement2D; i++)
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

} // namespace
