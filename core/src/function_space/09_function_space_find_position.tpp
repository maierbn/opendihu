#include "function_space/09_function_space_find_position.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"
#include "mesh/face_t.h"

namespace FunctionSpace
{

template<typename MeshType, typename BasisFunctionType>
bool FunctionSpaceFindPosition<MeshType,BasisFunctionType,Mesh::isStructured<MeshType>>::
findPosition(Vec3 point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,MeshType::dim()> &xi)
{
  const element_no_t nElements = this->nElementsLocal();
 
  // set starting no to 0 if it was not given and is thus arbitrarily initialized
  if (elementNo < 0 || elementNo >= nElements)
    elementNo = 0;
  
  // check if point is already in current element
  FunctionSpaceFindPosition<MeshType,BasisFunctionType> *functionSpace = this;

  if (ghostMeshNo != -1)
    functionSpace = ghostMesh_[ghostMeshNo].get();

  if (functionSpace->pointIsInElement(point, elementNo, xi))
  {
    return true;
  }
  
  // point is not in current element, consider the neighbouring elements and ghost meshes
  // neighbouring elements
  std::vector<std::pair<element_no_t,int>> neighbouringElements;   // (element no, ghostMeshNo)

  // set the neighbouring element nos, also considering ghost meshes
  getNeighbouringElements(elementNo, ghostMeshNo, neighbouringElements);

  for (int i = 0; i < neighbouringElements.size(); i++)
  {
    // leave out invalid neighbours (e.g. for element that are already at the border)
    if (neighbouringElements[i].first == -1)
      continue;

    element_no_t currentElementNo = neighbouringElements[i].first;

    // select ghost mesh or normal mesh
    FunctionSpaceFindPosition<MeshType,BasisFunctionType> *functionSpace = this;

    if (neighbouringElements[i].second != -1)
    {
      functionSpace = ghostMesh_[neighbouringElements[i].second].get();
    }

    // check if point is in current element
    if (functionSpace->pointIsInElement(point, currentElementNo, xi))
    {
      elementNo = currentElementNo;
      ghostMeshNo = neighbouringElements[i].second;
      return true;
    }
  }

  // search among all elements

  // look in every element, starting at elementNo-2
  element_no_t elementNoStart = (elementNo - 2 + nElements) % nElements;
  element_no_t elementNoEnd = (elementNo - 3 + nElements) % nElements;

  VLOG(3) << "elementNoStart: " << elementNoStart << ", elementNoEnd: " << elementNoEnd << ", nElements: " << nElements;

  LOG(WARNING) << "Could not find element that contains point " << point << " in neighbourhood of element " << elementNo << ", now testing all elements (no ghost elements).";
  for (element_no_t currentElementNo = elementNoStart; currentElementNo != elementNoEnd; currentElementNo++)
  {
    // restart with element 0
    if (currentElementNo == nElements)
    {
      currentElementNo = 0;
      if (elementNoEnd == currentElementNo)
        break;
    }

    if (this->pointIsInElement(point, currentElementNo, xi))
    {
      elementNo = currentElementNo;
      ghostMeshNo = -1;   // not a ghost mesh
      return true;
    }
  }
  return false;
}


template<typename MeshType, typename BasisFunctionType>
void FunctionSpaceFindPosition<MeshType,BasisFunctionType,Mesh::isStructured<MeshType>>::
setGhostMesh(Mesh::face_t face, std::shared_ptr<FunctionSpace<MeshType,BasisFunctionType>> ghostMesh)
{
  assert(0 <= face);
  assert(face < 6);
  ghostMesh_[face] == ghostMesh;
}

template<typename MeshType, typename BasisFunctionType>
void FunctionSpaceFindPosition<MeshType,BasisFunctionType,Mesh::isStructured<MeshType>>::
getNeighbouringElements(element_no_t elementNo, int ghostMeshNo, std::vector<std::pair<element_no_t,int>> &neighbouringElements)
{
  if (MeshType::dim() == 3)
  {
    neighbouringElements.resize(26);

    // get local coordinates of element
    std::array<int,3> coordinatesLocal;
    if (ghostMeshNo == -1)
    {
       coordinatesLocal = this->meshPartition_->getElementCoordinatesLocal(elementNo);
    }
    else if (ghostMeshNo == (int)Mesh::face_t::face0Minus)
    {
      assert(this->ghostMesh_[ghostMeshNo]);
      coordinatesLocal = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementCoordinatesLocal(elementNo);
      coordinatesLocal[0] = -1;
    }
    else if (ghostMeshNo == (int)Mesh::face_t::face0Plus)
    {
      assert(this->ghostMesh_[ghostMeshNo]);
      coordinatesLocal = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementCoordinatesLocal(elementNo);
      coordinatesLocal[0] = this->meshPartition_->nElementsLocal(0);
    }
    else if (ghostMeshNo == (int)Mesh::face_t::face1Minus)
    {
      assert(this->ghostMesh_[ghostMeshNo]);
      coordinatesLocal = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementCoordinatesLocal(elementNo);
      coordinatesLocal[1] = -1;
    }
    else if (ghostMeshNo == (int)Mesh::face_t::face1Plus)
    {
      assert(this->ghostMesh_[ghostMeshNo]);
      coordinatesLocal = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementCoordinatesLocal(elementNo);
      coordinatesLocal[1] = this->meshPartition_->nElementsLocal(1);
    }
    else if (ghostMeshNo == (int)Mesh::face_t::face2Minus)
    {
      assert(this->ghostMesh_[ghostMeshNo]);
      coordinatesLocal = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementCoordinatesLocal(elementNo);
      coordinatesLocal[2] = -1;
    }
    else if (ghostMeshNo == (int)Mesh::face_t::face2Plus)
    {
      assert(this->ghostMesh_[ghostMeshNo]);
      coordinatesLocal = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementCoordinatesLocal(elementNo);
      coordinatesLocal[2] = this->meshPartition_->nElementsLocal(2);
    }

    for (int x = -1; x <= +1; x++)
    {
      for (int y = -1; y <= +1; y++)
      {
        for (int z = -1; z <= +1; z++)
        {
          if (x == 0 && y == 0 && z == 0)
          {
            continue;
          }

          // compute index of neighbouringElements for the current element
          int index = (z+1)*9 + (y+1)*3 + (x+1);
          if (index >= 13)
            index--;

          neighbouringElements[index].first = -1;   // set to invalid

          coordinatesLocal += std::array<int,3>({x,y,z});

          if (coordinatesLocal[0] == -1)  // if at left boundary
          {
            if (ghostMesh_[(int)Mesh::face_t::face0Minus] != nullptr)
            {
              neighbouringElements[index].second = (int)Mesh::face_t::face0Minus;
              std::array<int,3> ghostMeshCoordinates({0,coordinatesLocal[1],coordinatesLocal[2]});
              neighbouringElements[index].first = this->ghostMesh_[(int)Mesh::face_t::face0Minus]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);
            }
            else
            {
              neighbouringElements[index].first = -1;   // set element no to invalid, because no ghost mesh was specified in that direction
            }
          }
          else if (coordinatesLocal[0] == this->meshPartition_->nElementsLocal(0))  // if at right boundary
          {
            if (ghostMesh_[(int)Mesh::face_t::face0Plus] != nullptr)
            {
              neighbouringElements[index].second = (int)Mesh::face_t::face0Plus;
              std::array<int,3> ghostMeshCoordinates({0,coordinatesLocal[1],coordinatesLocal[2]});
              neighbouringElements[index].first = this->ghostMesh_[(int)Mesh::face_t::face0Plus]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);
            }
            else
            {
              neighbouringElements[index].first = -1;   // set element no to invalid, because no ghost mesh was specified in that direction
            }
          }
          else if (coordinatesLocal[1] == 0)  // if at front boundary
          {
            if (ghostMesh_[(int)Mesh::face_t::face1Minus] != nullptr)
            {
              neighbouringElements[index].second = (int)Mesh::face_t::face1Minus;
              std::array<int,3> ghostMeshCoordinates({coordinatesLocal[0],0,coordinatesLocal[2]});
              neighbouringElements[index].first = this->ghostMesh_[(int)Mesh::face_t::face1Minus]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);
            }
            else
            {
              neighbouringElements[index].first = -1;   // set element no to invalid, because no ghost mesh was specified in that direction
            }
          }
          else if (coordinatesLocal[1] == this->meshPartition_->nElementsLocal(1))  // if at back boundary
          {
            if (ghostMesh_[(int)Mesh::face_t::face1Plus] != nullptr)
            {
              neighbouringElements[index].second = (int)Mesh::face_t::face1Plus;
              std::array<int,3> ghostMeshCoordinates({coordinatesLocal[0],0,coordinatesLocal[2]});
              neighbouringElements[index].first = this->ghostMesh_[(int)Mesh::face_t::face1Plus]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);
            }
            else
            {
              neighbouringElements[index].first = -1;   // set element no to invalid, because no ghost mesh was specified in that direction
            }
          }
          else if (coordinatesLocal[2] == 0)  // if at bottom boundary
          {
            if (ghostMesh_[(int)Mesh::face_t::face2Minus] != nullptr)
            {
              neighbouringElements[index].second = (int)Mesh::face_t::face2Minus;
              std::array<int,3> ghostMeshCoordinates({coordinatesLocal[0],coordinatesLocal[1],0});
              neighbouringElements[index].first = this->ghostMesh_[(int)Mesh::face_t::face2Minus]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);
            }
            else
            {
              neighbouringElements[index].first = -1;   // set element no to invalid, because no ghost mesh was specified in that direction
            }
          }
          else if (coordinatesLocal[2] == this->meshPartition_->nElementsLocal(2))  // if at top boundary
          {
            if (ghostMesh_[(int)Mesh::face_t::face2Plus] != nullptr)
            {
              neighbouringElements[index].second = (int)Mesh::face_t::face2Plus;
              std::array<int,3> ghostMeshCoordinates({coordinatesLocal[0],coordinatesLocal[1],0});
              neighbouringElements[index].first = this->ghostMesh_[(int)Mesh::face_t::face2Plus]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);
            }
            else
            {
              neighbouringElements[index].first = -1;   // set element no to invalid, because no ghost mesh was specified in that direction
            }
          }
          else if (0 <= coordinatesLocal[0] && coordinatesLocal[0] <= this->meshPartition_->nElementsLocal(0)
            && 0 <= coordinatesLocal[1] && coordinatesLocal[1] <= this->meshPartition_->nElementsLocal(1)
            && 0 <= coordinatesLocal[2] && coordinatesLocal[2] <= this->meshPartition_->nElementsLocal(2))
          {
            neighbouringElements[index].second = -1; // own mesh, no ghost mesh
            neighbouringElements[index].first = this->meshPartition_->getElementNoLocal(coordinatesLocal);
          }
        }
      }
    }
  }
  else assert(false);  // not implemented for D != 3
}

// unstructured
template<int D,typename BasisFunctionType>
bool FunctionSpaceFindPosition<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType,Mesh::UnstructuredDeformableOfDimension<D>>::
findPosition(Vec3 point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,D> &xi)
{
  const element_no_t nElements = this->nElementsLocal();

  // set starting no to 0 if it was not given and is thus arbitrarily initialized
  if (elementNo < 0 || elementNo >= nElements)
    elementNo = 0;

  // check if point is already in current element
  if (this->pointIsInElement(point, elementNo, xi))
  {
    return true;
  }

  // look in every element, starting at elementNo-2
  element_no_t elementNoStart = (elementNo - 2 + nElements) % nElements;
  element_no_t elementNoEnd = (elementNo - 3 + nElements) % nElements;

  VLOG(3) << "elementNoStart: " << elementNoStart << ", elementNoEnd: " << elementNoEnd << ", nElements: " << nElements;

  // search among all elements
  for (element_no_t currentElementNo = elementNoStart; currentElementNo != elementNoEnd; currentElementNo++)
  {
    // restart with element 0
    if (currentElementNo == nElements)
    {
      currentElementNo = 0;
      if (elementNoEnd == currentElementNo)
        break;
    }

    if (this->pointIsInElement(point, currentElementNo, xi))
    {
      elementNo = currentElementNo;
      ghostMeshNo = -1;   // not a ghost mesh
      return true;
    }
  }
  return false;
}


};  // namespace
