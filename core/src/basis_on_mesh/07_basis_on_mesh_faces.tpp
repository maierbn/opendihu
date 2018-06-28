#include "basis_on_mesh/07_basis_on_mesh_faces.h"

#include <Python.h>  // has to be the first included header
#include <array>
#include "easylogging++.h"

namespace BasisOnMesh
{

//! get all dof indices of a face for 1D
template<typename MeshType,typename BasisFunctionType>
void BasisOnMeshFaces<MeshType,BasisFunctionType,Mesh::isDim<1,MeshType>> ::
getFaceDofs(Mesh::face_t face, std::array<dof_no_t,BasisOnMeshBaseDim<0,BasisFunctionType>::nDofsPerElement()> &dofIndices)
{
  assert(face <= Mesh::face0Plus);
  switch(face)
  {
  case Mesh::face0Minus:
    dofIndices[0] = 0;
    break;

  case Mesh::face0Plus:
    dofIndices[0] = BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerElement()-1;
    break;
  }
}

//! get all dof indices of a face for 2D
template<typename MeshType,typename BasisFunctionType>
void BasisOnMeshFaces<MeshType,BasisFunctionType,Mesh::isDim<2,MeshType>> ::
getFaceDofs(Mesh::face_t face, std::array<dof_no_t,BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerElement()> &dofIndices)
{
  assert(face <= Mesh::face1Plus);

  //   1+
  //   _
  //0-|_|0+
  //   1-
  //
  dof_no_t currentDofIndex = 0;
  switch(face)
  {
  case Mesh::face0Minus:

    currentDofIndex = 0;
    for (int i=0; i<BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerElement(); i++)
    {
      dofIndices[i] = currentDofIndex;
      currentDofIndex += BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerElement();
    }
    break;

  case Mesh::face0Plus:

    currentDofIndex = BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerElement() - 1;
    for (int i=0; i<BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerElement(); i++)
    {
      dofIndices[i] = currentDofIndex;
      currentDofIndex += BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerElement();
    }
    break;

  case Mesh::face1Minus:

    currentDofIndex = 0;
    for (int i=0; i<BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerElement(); i++)
    {
      dofIndices[i] = currentDofIndex;
      currentDofIndex++;
    }
    break;

  case Mesh::face1Plus:

    currentDofIndex = BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerElement() - BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerElement();
    for (int i=0; i<BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerElement(); i++)
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
void BasisOnMeshFaces<MeshType,BasisFunctionType,Mesh::isDim<3,MeshType>> ::
getFaceDofs(Mesh::face_t face, std::array<dof_no_t,BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerElement()> &dofIndices)
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

  switch(face)
  {
  case Mesh::face0Minus:

    currentDofIndex = 0;
    for (int i=0; i<BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerElement(); i++)
    {
      dofIndices[i] = currentDofIndex;
      currentDofIndex += BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerElement();
    }
    break;

  case Mesh::face0Plus:

    currentDofIndex = BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerElement() - 1;
    for (int i=0; i<BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerElement(); i++)
    {
      dofIndices[i] = currentDofIndex;
      currentDofIndex += BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerElement();
    }
    break;

  case Mesh::face1Minus:

    currentDofIndex = 0;
    rowStartIndex = 0;
    index = 0;
    for (int i=0; i<BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerElement(); i++)
    {
      currentDofIndex = rowStartIndex;
      for (int j=0; j<BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerElement(); j++)
      {
        dofIndices[index++] = currentDofIndex;
        currentDofIndex++;
      }
      rowStartIndex += BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerElement();
    }
    break;

  case Mesh::face1Plus:

    currentDofIndex = 0;
    rowStartIndex = BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerElement() - BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerElement();
    index = 0;
    for (int i=0; i<BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerElement(); i++)
    {
      currentDofIndex = rowStartIndex;
      for (int j=0; j<BasisOnMeshBaseDim<1,BasisFunctionType>::nDofsPerElement(); j++)
      {
        dofIndices[index++] = currentDofIndex;
        currentDofIndex++;
      }
      rowStartIndex += BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerElement();
    }
    break;

  case Mesh::face2Minus:

    currentDofIndex = 0;
    for (int i=0; i<BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerElement(); i++)
    {
      dofIndices[i] = currentDofIndex;
      currentDofIndex++;
    }
    break;

  case Mesh::face2Plus:

    currentDofIndex = BasisOnMeshBaseDim<3,BasisFunctionType>::nDofsPerElement() - BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerElement();
    for (int i=0; i<BasisOnMeshBaseDim<2,BasisFunctionType>::nDofsPerElement(); i++)
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

};  // namespace