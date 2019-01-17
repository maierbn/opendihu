#include "function_space/07_function_space_faces.h"

#include <Python.h>  // has to be the first included header
#include <array>
#include "easylogging++.h"

namespace FunctionSpace
{

//! get all dof indices of a face for 1D
template<typename MeshType,typename BasisFunctionType>
void FunctionSpaceFaces<MeshType,BasisFunctionType,Mesh::isDim<1,MeshType>> ::
getFaceDofs(Mesh::face_t face, std::array<dof_no_t,FunctionSpaceBaseDim<0,BasisFunctionType>::nDofsPerElement()> &dofIndices)
{
  assert(face <= Mesh::face0Plus);
  switch(face)
  {
  case Mesh::face0Minus:
    dofIndices[0] = 0;
    break;

  case Mesh::face0Plus:
    dofIndices[0] = FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerElement()-1;
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
  switch(face)
  {
  case Mesh::face0Minus:

    currentDofIndex = 0;
    for (int i=0; i<FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerElement(); i++)
    {
      dofIndices[i] = currentDofIndex;
      currentDofIndex += FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerElement();
    }
    break;

  case Mesh::face0Plus:

    currentDofIndex = FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerElement() - 1;
    for (int i=0; i<FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerElement(); i++)
    {
      dofIndices[i] = currentDofIndex;
      currentDofIndex += FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerElement();
    }
    break;

  case Mesh::face1Minus:

    currentDofIndex = 0;
    for (int i=0; i<FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerElement(); i++)
    {
      dofIndices[i] = currentDofIndex;
      currentDofIndex++;
    }
    break;

  case Mesh::face1Plus:

    currentDofIndex = FunctionSpaceBaseDim<2,BasisFunctionType>::nDofsPerElement() - FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerElement();
    for (int i=0; i<FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerElement(); i++)
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

  switch(face)
  {
  case Mesh::face0Minus:

    currentDofIndex = 0;
    for (int i=0; i<FunctionSpaceBaseDim<2,BasisFunctionType>::nDofsPerElement(); i++)
    {
      dofIndices[i] = currentDofIndex;
      currentDofIndex += FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerElement();
    }
    break;

  case Mesh::face0Plus:

    currentDofIndex = FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerElement() - 1;
    for (int i=0; i<FunctionSpaceBaseDim<2,BasisFunctionType>::nDofsPerElement(); i++)
    {
      dofIndices[i] = currentDofIndex;
      currentDofIndex += FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerElement();
    }
    break;

  case Mesh::face1Minus:

    currentDofIndex = 0;
    rowStartIndex = 0;
    index = 0;
    for (int i=0; i<FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerElement(); i++)
    {
      currentDofIndex = rowStartIndex;
      for (int j=0; j<FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerElement(); j++)
      {
        dofIndices[index++] = currentDofIndex;
        currentDofIndex++;
      }
      rowStartIndex += FunctionSpaceBaseDim<2,BasisFunctionType>::nDofsPerElement();
    }
    break;

  case Mesh::face1Plus:

    currentDofIndex = 0;
    rowStartIndex = FunctionSpaceBaseDim<2,BasisFunctionType>::nDofsPerElement() - FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerElement();
    index = 0;
    for (int i=0; i<FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerElement(); i++)
    {
      currentDofIndex = rowStartIndex;
      for (int j=0; j<FunctionSpaceBaseDim<1,BasisFunctionType>::nDofsPerElement(); j++)
      {
        dofIndices[index++] = currentDofIndex;
        currentDofIndex++;
      }
      rowStartIndex += FunctionSpaceBaseDim<2,BasisFunctionType>::nDofsPerElement();
    }
    break;

  case Mesh::face2Minus:

    currentDofIndex = 0;
    for (int i=0; i<FunctionSpaceBaseDim<2,BasisFunctionType>::nDofsPerElement(); i++)
    {
      dofIndices[i] = currentDofIndex;
      currentDofIndex++;
    }
    break;

  case Mesh::face2Plus:

    currentDofIndex = FunctionSpaceBaseDim<3,BasisFunctionType>::nDofsPerElement() - FunctionSpaceBaseDim<2,BasisFunctionType>::nDofsPerElement();
    for (int i=0; i<FunctionSpaceBaseDim<2,BasisFunctionType>::nDofsPerElement(); i++)
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
