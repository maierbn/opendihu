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
    for (int nodalDofIndex = 0; nodalDofIndex != nDofsPerNode; nodalDofIndex++)
    {
      dofIndices[nodalDofIndex] = nodalDofIndex;
    }
    break;

  case Mesh::face0Plus:
    for (int nodalDofIndex = 0; nodalDofIndex != nDofsPerNode; nodalDofIndex++)
    {
      dofIndices[nodalDofIndex] = nDofsPerElement-nDofsPerNode + nodalDofIndex;
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
getFaceDofs(Mesh::face_t face, std::array<dof_no_t,FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerElement()> &dofIndices)
{
  assert(face <= Mesh::face1Plus);

  //   1+
  //   _
  //0-|_|0+
  //   1-
  //
  dof_no_t currentDofIndex = 0;
  const int nDofsPerElement = FunctionSpaceBaseDim<MeshType::dim(),BasisFunctionType>::nDofsPerElement();
  const int nDofsPerNode = FunctionSpaceBaseDim<MeshType::dim(),BasisFunctionType>::nDofsPerNode();
  //const int nDofsPerNode1D = FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerNode();
  switch(face)
  {
  case Mesh::face0Minus:

    currentDofIndex = 0;
    for (int i = 0; i < nDofsPerElement; i++)
    {
      for (int nodalDofIndex = 0; nodalDofIndex != nDofsPerNode; nodalDofIndex++, i++)
      {
        dofIndices[i] = currentDofIndex + nodalDofIndex;
      }
      currentDofIndex += nDofsPerElement;
    }
    break;

  case Mesh::face0Plus:

    currentDofIndex = nDofsPerElement - nDofsPerNode;
    for (int i = 0; i < nDofsPerElement; i++)
    {
      for (int nodalDofIndex = 0; nodalDofIndex != nDofsPerNode; nodalDofIndex++, i++)
      {
        dofIndices[i] = currentDofIndex + nodalDofIndex;
      }
      currentDofIndex += nDofsPerElement;
    }
    break;

  case Mesh::face1Minus:

    currentDofIndex = 0;
    for (int i = 0; i < nDofsPerElement; i++)
    {
      dofIndices[i] = currentDofIndex;
      currentDofIndex++;
    }
    break;

  case Mesh::face1Plus:

    currentDofIndex = FunctionSpaceBaseDim<2,BasisFunctionType>::nDofsPerElement() - nDofsPerElement;
    for (int i = 0; i < nDofsPerElement; i++)
    {
      dofIndices[i] = currentDofIndex;
      currentDofIndex++;
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
getFaceDofs(Mesh::face_t face, std::array<dof_no_t,FunctionSpaceBaseDim<2,BasisFunctionType>::nDofsPerElement()> &dofIndices)
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

  const int nDofsPerElement2D = FunctionSpaceBaseDim<2,BasisFunctionType>::nDofsPerElement();
  const int nDofsPerElement1D = FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerElement();
  const int nDofsPerNode = FunctionSpaceBaseDim<MeshType::dim(),BasisFunctionType>::nDofsPerNode();

  switch(face)
  {
  case Mesh::face0Minus:

    currentDofIndex = 0;
    for (int i = 0; i < nDofsPerElement2D; i++)
    {
      for (int nodalDofIndex = 0; nodalDofIndex != nDofsPerNode; nodalDofIndex++, i++)
      {
        dofIndices[i] = currentDofIndex + nodalDofIndex;
      }
      currentDofIndex += nDofsPerElement1D;
    }
    break;

  case Mesh::face0Plus:

    currentDofIndex = nDofsPerElement1D - nDofsPerNode;
    for (int i = 0; i < nDofsPerElement2D; i++)
    {
      for (int nodalDofIndex = 0; nodalDofIndex != nDofsPerNode; nodalDofIndex++, i++)
      {
        dofIndices[i] = currentDofIndex + nodalDofIndex;
      }
      currentDofIndex += nDofsPerElement1D;
    }
    break;

  case Mesh::face1Minus:

    currentDofIndex = 0;
    rowStartIndex = 0;
    index = 0;
    for (int i = 0; i < nDofsPerElement1D; i++)
    {
      currentDofIndex = rowStartIndex;
      for (int j = 0; j < nDofsPerElement1D; j++)
      {
        dofIndices[index++] = currentDofIndex;
        currentDofIndex++;
      }
      rowStartIndex += nDofsPerElement2D;
    }
    break;

  case Mesh::face1Plus:

    currentDofIndex = 0;
    rowStartIndex = nDofsPerElement2D - nDofsPerElement1D;
    index = 0;
    for (int i = 0; i < nDofsPerElement1D; i++)
    {
      currentDofIndex = rowStartIndex;
      for (int j = 0; j < nDofsPerElement1D; j++)
      {
        dofIndices[index++] = currentDofIndex;
        currentDofIndex++;
      }
      rowStartIndex += nDofsPerElement2D;
    }
    break;

  case Mesh::face2Minus:

    currentDofIndex = 0;
    for (int i = 0; i < nDofsPerElement2D; i++)
    {
      dofIndices[i] = currentDofIndex;
      currentDofIndex++;
    }
    break;

  case Mesh::face2Plus:

    currentDofIndex = FunctionSpaceBaseDim<3,BasisFunctionType>::nDofsPerElement() - nDofsPerElement2D;
    for (int i = 0; i < nDofsPerElement2D; i++)
    {
      dofIndices[i] = currentDofIndex;
      currentDofIndex++;
    }
    break;

  default:
    assert(false);
    break;
  }
}

} // namespace
