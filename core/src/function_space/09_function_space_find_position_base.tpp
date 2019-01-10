#include "function_space/09_function_space_find_position_base.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"
#include "mesh/face_t.h"

namespace FunctionSpace
{

// structured mesh
template<typename MeshType, typename BasisFunctionType>
bool FunctionSpaceStructuredFindPositionBase<MeshType,BasisFunctionType>::
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
    FunctionSpaceStructuredFindPositionBase<MeshType,BasisFunctionType> *functionSpace = this;

    if (ghostMeshNo != -1)
      functionSpace = ghostMesh_[ghostMeshNo].get();

    if (functionSpace->pointIsInElement(point, elementNo, xi))
    {

      // debugging output
      if (VLOG_IS_ON(2) && MeshType::dim() == 3)
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
    if (this->checkNeighbouringElements(point, elementNo, ghostMeshNo, xi))
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
    else
    {
      VLOG(4) << " no (xi=" << xi << ")";
    }
  }

  // if point was still not found, search in ghost meshes
  for (int face = (int)Mesh::face_t::face0Minus; face <= (int)Mesh::face_t::face2Plus; face++)
  {
    VLOG(3) << "consider ghost mesh " << Mesh::getString((Mesh::face_t)face);
    if (ghostMesh_[face] != nullptr)
    {
      VLOG(3) << "   ghost mesh " << Mesh::getString((Mesh::face_t)face) << " is set";
      if (ghostMesh_[face]->findPosition(point, elementNo, ghostMeshNo, xi, false))
      {
        VLOG(3) << "   point found in ghost mesh " << Mesh::getString((Mesh::face_t)face) << ", element " << elementNo << ", xi " << xi;
        ghostMeshNo = face;
        return true;
      }
      else
      {
        VLOG(3) << "   not found";
      }
    }
    else
    {
      VLOG(3) << "   ghost mesh " << Mesh::getString((Mesh::face_t)face) << " is not set";
    }
  }

  VLOG(1) << "Could not find any containing element (streamline ends)";
  return false;
}

template<typename MeshType, typename BasisFunctionType>
void FunctionSpaceStructuredFindPositionBase<MeshType,BasisFunctionType>::
setGhostMesh(Mesh::face_t face, const std::shared_ptr<FunctionSpace<MeshType,BasisFunctionType>> ghostMesh)
{
  assert(0 <= face);
  assert(face < 6);
  ghostMesh_[face] = ghostMesh;
  VLOG(1) << "set ghost mesh for face " << Mesh::getString((Mesh::face_t)face) << " to " << (ghostMesh == nullptr? " null" : "x");
}

template<typename MeshType, typename BasisFunctionType>
void FunctionSpaceStructuredFindPositionBase<MeshType,BasisFunctionType>::
debugOutputGhostMeshSet()
{
  VLOG(1) << "ghost mesh 0- is " << (ghostMesh_[(int)Mesh::face_t::face0Minus] == nullptr? "not" : "") << " set";
  VLOG(1) << "ghost mesh 0+ is " << (ghostMesh_[(int)Mesh::face_t::face0Plus] == nullptr? "not" : "") << " set";
  VLOG(1) << "ghost mesh 1- is " << (ghostMesh_[(int)Mesh::face_t::face1Minus] == nullptr? "not" : "") << " set";
  VLOG(1) << "ghost mesh 1+ is " << (ghostMesh_[(int)Mesh::face_t::face1Plus] == nullptr? "not" : "") << " set";
  VLOG(1) << "ghost mesh 2- is " << (ghostMesh_[(int)Mesh::face_t::face2Minus] == nullptr? "not" : "") << " set";
  VLOG(1) << "ghost mesh 2+ is " << (ghostMesh_[(int)Mesh::face_t::face2Plus] == nullptr? "not" : "") << " set";
}

};  // namespace
