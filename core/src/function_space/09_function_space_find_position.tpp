#include "function_space/09_function_space_find_position.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"
#include "mesh/face_t.h"

namespace FunctionSpace
{

// structured mesh
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
    VLOG(2) << "findPosition: startSearchInCurrentElement";

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

    VLOG(2) << "point is not in current element, now check neighbouring elements";

    // set the neighbouring element nos, also considering ghost meshes
    if (checkNeighbouringElements(point, elementNo, ghostMeshNo, xi))
    {
      return true;
    }
    else
    {
      VLOG(2) << "point was also not found among neighbouring elements (including ghost meshes), xi: " << xi << ", now check all elements";
    }
  }

  // search among all elements

  // look in every element, starting at elementNo-2
  element_no_t elementNoStart = (elementNo - 2 + nElements) % nElements;
  element_no_t elementNoEnd = (elementNo - 3 + nElements) % nElements;

  VLOG(3) << "elementNoStart: " << elementNoStart << ", elementNoEnd: " << elementNoEnd << ", nElements: " << nElements
    << "(" << this->meshPartition_->nElementsLocal(0) << "x" << this->meshPartition_->nElementsLocal(1) << "x" << this->meshPartition_->nElementsLocal(2) << ")";

  for (element_no_t currentElementNo = elementNoStart; currentElementNo != elementNoEnd; currentElementNo++)
  {
    // restart with element 0
    if (currentElementNo == nElements)
    {
      currentElementNo = 0;
      if (elementNoEnd == currentElementNo)
        break;
    }

    VLOG(4) << "check element " << currentElementNo;

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
setGhostMesh(Mesh::face_t face, const std::shared_ptr<FunctionSpace<MeshType,BasisFunctionType>> ghostMesh)
{
  assert(0 <= face);
  assert(face < 6);
  ghostMesh_[face] = ghostMesh;
  VLOG(1) << "set ghost mesh for face " << Mesh::getString((Mesh::face_t)face) << " to " << (ghostMesh == nullptr? " null" : "x");
}

template<typename MeshType, typename BasisFunctionType>
void FunctionSpaceFindPosition<MeshType,BasisFunctionType,Mesh::isStructured<MeshType>>::
debugOutputGhostMeshSet()
{
  VLOG(1) << "ghost mesh 0- is " << (ghostMesh_[(int)Mesh::face_t::face0Minus] == nullptr? "not" : "") << " set";
  VLOG(1) << "ghost mesh 0+ is " << (ghostMesh_[(int)Mesh::face_t::face0Plus] == nullptr? "not" : "") << " set";
  VLOG(1) << "ghost mesh 1- is " << (ghostMesh_[(int)Mesh::face_t::face1Minus] == nullptr? "not" : "") << " set";
  VLOG(1) << "ghost mesh 1+ is " << (ghostMesh_[(int)Mesh::face_t::face1Plus] == nullptr? "not" : "") << " set";
  VLOG(1) << "ghost mesh 2- is " << (ghostMesh_[(int)Mesh::face_t::face2Minus] == nullptr? "not" : "") << " set";
  VLOG(1) << "ghost mesh 2+ is " << (ghostMesh_[(int)Mesh::face_t::face2Plus] == nullptr? "not" : "") << " set";
}

template<typename MeshType, typename BasisFunctionType>
bool FunctionSpaceFindPosition<MeshType,BasisFunctionType,Mesh::isStructured<MeshType>>::
checkNeighbouringElements(const Vec3 &point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,MeshType::dim()> &xi)
{
  VLOG(1) << "getNeighbouringElements(elementNo = " << elementNo << ", ghostMeshNo = " << ghostMeshNo << ", initial xi = " << xi;

  if (MeshType::dim() == 3)
  {
    static std::array<int,3> xOffset;
    static std::array<int,3> yOffset;
    static std::array<int,3> zOffset;

    // x direction
    if (xi[0] < 0)
    {
      xOffset = {-1, 0, 1};
    }
    else if (xi[0] > 1.0)
    {
      xOffset = {1, 0, -1};
    }
    else if (xi[0] > 0.5)
    {
      xOffset = {0, 1, -1};
    }
    else
    {
      xOffset = {0, -1, 1};
    }

    // y direction
    if (xi[1] < 0)
    {
      yOffset = {-1, 0, 1};
    }
    else if (xi[1] > 1.0)
    {
      yOffset = {1, 0, -1};
    }
    else if (xi[1] > 0.5)
    {
      yOffset = {0, 1, -1};
    }
    else
    {
      yOffset = {0, -1, 1};
    }

    // z direction
    if (xi[2] < 0)
    {
      zOffset = {-1, 0, 1};
    }
    else if (xi[2] > 1.0)
    {
      zOffset = {1, 0, -1};
    }
    else if (xi[2] > 0.5)
    {
      zOffset = {0, 1, -1};
    }
    else
    {
      zOffset = {0, -1, 1};
    }

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
      << "," << this->meshPartition_->nElementsLocal(2) << "] coordinatesLocal: " << coordinatesLocal
      << " interation z in " << zOffset << ", y in " << yOffset << " x in " << xOffset << "";

    //debugOutputGhostMeshSet();

    std::array<int,3> neighbourCoordinatesLocal;
    element_no_t neighbourElementNo;

    FunctionSpaceFindPosition<MeshType,BasisFunctionType> *functionSpace = this;

    for (int zIndex = 0; zIndex != 3; zIndex++)
    {
      int z = zOffset[zIndex];
      for (int yIndex = 0; yIndex != 3; yIndex++)
      {
        int y = yOffset[yIndex];
        for (int xIndex = 0; xIndex != 3; xIndex++)
        {
          int x = xOffset[xIndex];

          // do not handle the center point
          if (x == 0 && y == 0 && z == 0)
          {
            continue;
          }

          neighbourElementNo = -1;   // set to invalid

          neighbourCoordinatesLocal = coordinatesLocal + std::array<int,3>({x,y,z});

          VLOG(1) << "(x,y,z) = (" << x << "," << y << "," << z << "), neighbourCoordinatesLocal: " << neighbourCoordinatesLocal;

          if (neighbourCoordinatesLocal[2] > this->meshPartition_->nElementsLocal(2) || neighbourCoordinatesLocal[2] < -1
            || neighbourCoordinatesLocal[1] > this->meshPartition_->nElementsLocal(1) || neighbourCoordinatesLocal[1] < -1
            || neighbourCoordinatesLocal[0] > this->meshPartition_->nElementsLocal(0) || neighbourCoordinatesLocal[0] < -1)
          {
            VLOG(1) << "outside ghost layer";
            continue;
          }

          // do not consider diagonal ghost meshes, e.g. x+/y+
          int nGhostTargets = 0;
          if (neighbourCoordinatesLocal[2] == this->meshPartition_->nElementsLocal(2) || neighbourCoordinatesLocal[2] == -1)
          {
            // if top and bottom ghost meshes were not provided, continue
            if (ghostMesh_[(int)Mesh::face_t::face2Minus] == nullptr && ghostMesh_[(int)Mesh::face_t::face2Plus] == nullptr)
            {
              VLOG(1) << "z+/z- ghost meshes not set, do not consider respective neighbours";
              continue;
            }
            nGhostTargets++;
          }
          if (neighbourCoordinatesLocal[1] == this->meshPartition_->nElementsLocal(1) || neighbourCoordinatesLocal[1] == -1)
          {
            nGhostTargets++;
          }
          if (neighbourCoordinatesLocal[0] == this->meshPartition_->nElementsLocal(0) || neighbourCoordinatesLocal[0] == -1)
          {
            nGhostTargets++;
          }
          if (nGhostTargets > 1)
          {
            VLOG(1) << "edge ghost element";
            continue;
          }

          functionSpace = this;

          if (neighbourCoordinatesLocal[0] == -1)  // if at left boundary
          {
            if (ghostMesh_[(int)Mesh::face_t::face0Minus] != nullptr)
            {
              ghostMeshNo = (int)Mesh::face_t::face0Minus;
              functionSpace = ghostMesh_[ghostMeshNo].get();
              std::array<int,3> ghostMeshCoordinates({0,neighbourCoordinatesLocal[1],neighbourCoordinatesLocal[2]});
              neighbourElementNo = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);

              VLOG(1) << "0- neighbourElementNo: " << neighbourElementNo << ", ghostMeshCoordinates: " << ghostMeshCoordinates
                << " / (" << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(2) << ")";
            }
            else
            {
              // at this position outside of the normal subdomain no ghost mesh was specified
              continue;
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
              VLOG(1) << "0+ neighbourElementNo: " << neighbourElementNo << ", ghostMeshCoordinates: " << ghostMeshCoordinates
                << " / (" << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(2) << ")";
            }
            else
            {
              // at this position outside of the normal subdomain no ghost mesh was specified
              continue;
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
              VLOG(1) << "1- neighbourElementNo: " << neighbourElementNo << ", ghostMeshCoordinates: " << ghostMeshCoordinates
                << " / (" << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(2) << ")";
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
              VLOG(1) << "1+ neighbourElementNo: " << neighbourElementNo << ", ghostMeshCoordinates: " << ghostMeshCoordinates
                << " / (" << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(2) << ")";
            }
            else
            {
              // at this position outside of the normal subdomain no ghost mesh was specified
              continue;
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
              VLOG(1) << "2- neighbourElementNo: " << neighbourElementNo << ", ghostMeshCoordinates: " << ghostMeshCoordinates
                << " / (" << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(2) << ")";
            }
            else
            {
              // at this position outside of the normal subdomain no ghost mesh was specified
              continue;
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
              VLOG(1) << "2+ neighbourElementNo: " << neighbourElementNo << ", ghostMeshCoordinates: " << ghostMeshCoordinates
                << " / (" << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(2) << ")";
            }
            else
            {
              // at this position outside of the normal subdomain no ghost mesh was specified
              continue;
            }
          }
          else if (0 <= neighbourCoordinatesLocal[0] && neighbourCoordinatesLocal[0] <= this->meshPartition_->nElementsLocal(0)
            && 0 <= neighbourCoordinatesLocal[1] && neighbourCoordinatesLocal[1] <= this->meshPartition_->nElementsLocal(1)
            && 0 <= neighbourCoordinatesLocal[2] && neighbourCoordinatesLocal[2] <= this->meshPartition_->nElementsLocal(2))
          {
            neighbourElementNo = this->meshPartition_->getElementNoLocal(neighbourCoordinatesLocal);
            functionSpace = this;
            ghostMeshNo = -1;
            VLOG(1) << "normal, neighbourElementNo: " << neighbourElementNo << " nElementsLocal: " << this->nElementsLocal() << "=" << functionSpace->nElementsLocal();
          }

          // check if point is in current neighbour element
          if (neighbourElementNo != -1)
          {
            if (functionSpace->pointIsInElement(point, neighbourElementNo, xi))
            {
              elementNo = neighbourElementNo;
              return true;
            }
            else
            {
              VLOG(1) << "  point is not in element " << neighbourElementNo << " (xi: " << xi << ").";
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
