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
findPosition(Vec3 point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,MeshType::dim()> &xi, bool startSearchInCurrentElement)
{
  const element_no_t nElements = this->nElementsLocal();
 
  // set starting no to 0 if it was not given and is thus arbitrarily initialized
  if (elementNo < 0 || elementNo >= nElements)
    elementNo = 0;

  if (startSearchInCurrentElement)
  {

    // check if point is already in current element
    FunctionSpaceFindPosition<MeshType,BasisFunctionType> *functionSpace = this;

    if (ghostMeshNo != -1)
      functionSpace = ghostMesh_[ghostMeshNo].get();

    if (functionSpace->pointIsInElement(point, elementNo, xi))
    {

      // debugging output
      if (VLOG_IS_ON(2))
      {
        // check for element size
        const int nDofsPerElement = FunctionSpace<MeshType,BasisFunctionType>::nDofsPerElement();
        std::array<Vec3,nDofsPerElement> elementalNodePositions;
        this->geometryField().getElementValues(elementNo, elementalNodePositions);

        // get bounding box of element
        double xMin = elementalNodePositions[0][0];
        double xMax = elementalNodePositions[0][0];
        double yMin = elementalNodePositions[0][1];
        double yMax = elementalNodePositions[0][1];
        double zMin = elementalNodePositions[0][2];
        double zMax = elementalNodePositions[0][2];

        for (Vec3 &elementalNodePosition : elementalNodePositions)
        {
          xMin = std::min(xMin, elementalNodePosition[0]);
          xMax = std::max(xMax, elementalNodePosition[0]);
          yMin = std::min(yMin, elementalNodePosition[1]);
          yMax = std::max(yMax, elementalNodePosition[1]);
          zMin = std::min(zMin, elementalNodePosition[2]);
          zMax = std::max(zMax, elementalNodePosition[2]);
        }

        //double xLength = xMax - xMin;
        //double yLength = yMax - yMin;
        //douuble zLength = zMax - zMin;

        VLOG(1) << "point " << point << " is in element " << elementNo << ", which has "
          << "bounding box x: [" << xMin << "," << xMax << "], y: [" << yMin << "," << yMax << "], z: [" << zMin << "," << zMax << "]";
      }

      return true;
    }

    // point is not in current element, consider the neighbouring elements and ghost meshes

    // set the neighbouring element nos, also considering ghost meshes
    if (checkNeighbouringElements(point, elementNo, ghostMeshNo, xi))
      return true;
  }

  // search among all elements

  // look in every element, starting at elementNo-2
  element_no_t elementNoStart = (elementNo - 2 + nElements) % nElements;
  element_no_t elementNoEnd = (elementNo - 3 + nElements) % nElements;

  VLOG(3) << "elementNoStart: " << elementNoStart << ", elementNoEnd: " << elementNoEnd << ", nElements: " << nElements;

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
      if (startSearchInCurrentElement)
      {
        LOG(WARNING) << "Could not find element that contains point " << point << " in neighbourhood of element " << elementNo
          << ", tested all elements (no ghost elements) and found element " << currentElementNo;
      }

      elementNo = currentElementNo;
      ghostMeshNo = -1;   // not a ghost mesh

      return true;
    }
  }

  VLOG(1) << "Could not find any containing element (streamline ends)";
  return false;
}


template<typename MeshType, typename BasisFunctionType>
void FunctionSpaceFindPosition<MeshType,BasisFunctionType,Mesh::isStructured<MeshType>>::
setGhostMesh(Mesh::face_t face, std::shared_ptr<FunctionSpace<MeshType,BasisFunctionType>> ghostMesh)
{
  assert(0 <= face);
  assert(face < 6);
  ghostMesh_[face] = ghostMesh;
}

template<typename MeshType, typename BasisFunctionType>
bool FunctionSpaceFindPosition<MeshType,BasisFunctionType,Mesh::isStructured<MeshType>>::
checkNeighbouringElements(const Vec3 &point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,MeshType::dim()> &xi)
{
  VLOG(1) << "getNeighbouringElements(elementNo = " << elementNo << ", ghostMeshNO = " << ghostMeshNo << ")";

  if (MeshType::dim() == 3)
  {
    // this is implemented in a generic way such that we could also set different sizes of the box of considered neighbours here (however, not tested)
    const int xLeft = -1;
    const int xRight = 1;
    const int yLeft = -1;
    const int yRight = 1;
    const int zLeft = -1;
    const int zRight = 1;

    // set loop boundaries, such that the last successful case is first
    int zBegin = -1;
    int zEnd = +2;
    int zIncrement = +1;
    if (targetZ_ == +1)
    {
      zBegin = +1;
      zEnd = -2;
      zIncrement = -1;
    }

    int yBegin = -1;
    int yEnd = +2;
    int yIncrement = +1;
    if (targetY_ == +1)
    {
      yBegin = +1;
      yEnd = -2;
      yIncrement = -1;
    }

    int xBegin = -1;
    int xEnd = +2;
    int xIncrement = +1;
    if (targetX_ == +1)
    {
      xBegin = +1;
      xEnd = -2;
      xIncrement = -1;
    }

    int nNeighbours = (xRight-xLeft+1) * (yRight-yLeft+1) * (zRight-zLeft+1) - 1;
    assert(nNeighbours == 26);
    //neighbouringElements.resize(nNeighbours);

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

    VLOG(1) << "nElementsLocal: [" << this->meshPartition_->nElementsLocal(0) << "," << this->meshPartition_->nElementsLocal(1)
      << "," << this->meshPartition_->nElementsLocal(2) << "] coordinatesLocal: " << coordinatesLocal;

    std::array<int,3> neighbourCoordinatesLocal;
    element_no_t neighbourElementNo;

    for (int z = zBegin; z != zEnd; z += zIncrement)
    {
      for (int y = yBegin; y != yEnd; y += yIncrement)
      {
        for (int x = xBegin; x != xEnd; x += xIncrement)
        {
          // do not handle the center point
          if (x == 0 && y == 0 && z == 0)
          {
            continue;
          }

          // compute index of neighbouringElements for the current element
          int index = (z-zLeft)*(xRight-xLeft+1)*(yRight-yLeft+1) + (y-yLeft)*(xRight-xLeft+1) + (x-xLeft);
          if (index >= nNeighbours/2)
            index--;

          neighbourElementNo = -1;   // set to invalid

          neighbourCoordinatesLocal = coordinatesLocal + std::array<int,3>({x,y,z});

          FunctionSpaceFindPosition<MeshType,BasisFunctionType> *functionSpace = this;

          VLOG(2) << "(x,y,z) = (" << x << "," << y << "," << z << "), index = " << index << ", neighbourCoordinatesLocal: " << neighbourCoordinatesLocal;

          if (neighbourCoordinatesLocal[0] == -1)  // if at left boundary
          {
            if (ghostMesh_[(int)Mesh::face_t::face0Minus] != nullptr)
            {
              ghostMeshNo = (int)Mesh::face_t::face0Minus;
              functionSpace = ghostMesh_[ghostMeshNo].get();
              std::array<int,3> ghostMeshCoordinates({0,neighbourCoordinatesLocal[1],neighbourCoordinatesLocal[2]});
              neighbourElementNo = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);
            }
          }
          else if (neighbourCoordinatesLocal[0] == this->meshPartition_->nElementsLocal(0))  // if at right boundary
          {
            if (ghostMesh_[(int)Mesh::face_t::face0Plus] != nullptr)
            {
              ghostMeshNo = (int)Mesh::face_t::face0Plus;
              functionSpace = ghostMesh_[ghostMeshNo].get();
              std::array<int,3> ghostMeshCoordinates({0,neighbourCoordinatesLocal[1],neighbourCoordinatesLocal[2]});
              neighbourElementNo = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);
            }
          }
          else if (neighbourCoordinatesLocal[1] == -1)  // if at front boundary
          {
            if (ghostMesh_[(int)Mesh::face_t::face1Minus] != nullptr)
            {
              ghostMeshNo = (int)Mesh::face_t::face1Minus;
              functionSpace = ghostMesh_[ghostMeshNo].get();
              std::array<int,3> ghostMeshCoordinates({neighbourCoordinatesLocal[0],0,neighbourCoordinatesLocal[2]});
              neighbourElementNo = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);
            }
          }
          else if (neighbourCoordinatesLocal[1] == this->meshPartition_->nElementsLocal(1))  // if at back boundary
          {
            if (ghostMesh_[(int)Mesh::face_t::face1Plus] != nullptr)
            {
              ghostMeshNo = (int)Mesh::face_t::face1Plus;
              functionSpace = ghostMesh_[ghostMeshNo].get();
              std::array<int,3> ghostMeshCoordinates({neighbourCoordinatesLocal[0],0,neighbourCoordinatesLocal[2]});
              neighbourElementNo = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);
            }
          }
          else if (neighbourCoordinatesLocal[2] == -1)  // if at bottom boundary
          {
            if (ghostMesh_[(int)Mesh::face_t::face2Minus] != nullptr)
            {
              ghostMeshNo = (int)Mesh::face_t::face2Minus;
              functionSpace = ghostMesh_[ghostMeshNo].get();
              std::array<int,3> ghostMeshCoordinates({neighbourCoordinatesLocal[0],neighbourCoordinatesLocal[1],0});
              neighbourElementNo = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);
            }
          }
          else if (neighbourCoordinatesLocal[2] == this->meshPartition_->nElementsLocal(2))  // if at top boundary
          {
            if (ghostMesh_[(int)Mesh::face_t::face2Plus] != nullptr)
            {
              ghostMeshNo = (int)Mesh::face_t::face2Plus;
              functionSpace = ghostMesh_[ghostMeshNo].get();
              std::array<int,3> ghostMeshCoordinates({neighbourCoordinatesLocal[0],neighbourCoordinatesLocal[1],0});
              neighbourElementNo = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);
            }
          }
          else if (0 <= neighbourCoordinatesLocal[0] && neighbourCoordinatesLocal[0] <= this->meshPartition_->nElementsLocal(0)
            && 0 <= neighbourCoordinatesLocal[1] && neighbourCoordinatesLocal[1] <= this->meshPartition_->nElementsLocal(1)
            && 0 <= neighbourCoordinatesLocal[2] && neighbourCoordinatesLocal[2] <= this->meshPartition_->nElementsLocal(2))
          {
            neighbourElementNo = this->meshPartition_->getElementNoLocal(neighbourCoordinatesLocal);
          }

          // check if point is in current neighbour element
          if (neighbourElementNo != -1)
          {
            if (functionSpace->pointIsInElement(point, neighbourElementNo, xi))
            {
              targetX_ = x;
              targetY_ = y;
              targetZ_ = z;
              elementNo = neighbourElementNo;
              return true;
            }
          }

        }  // z
      }  // y
    }   // x
  }  // dim == 3
  else assert(false);  // not implemented for D != 3
  return false;
}

// unstructured
template<int D,typename BasisFunctionType>
bool FunctionSpaceFindPosition<Mesh::UnstructuredDeformableOfDimension<D>,BasisFunctionType,Mesh::UnstructuredDeformableOfDimension<D>>::
findPosition(Vec3 point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,D> &xi, bool startSearchInCurrentElement)
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
