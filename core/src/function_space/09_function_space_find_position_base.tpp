#include "function_space/09_function_space_find_position_base.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"
#include "mesh/face_t.h"
#include "control/dihu_context.h"

namespace FunctionSpace
{

// structured mesh
template<typename MeshType, typename BasisFunctionType>
bool FunctionSpaceStructuredFindPositionBase<MeshType,BasisFunctionType>::
findPosition(Vec3 point, element_no_t &elementNoLocal, int &ghostMeshNo, std::array<double,MeshType::dim()> &xi, bool startSearchInCurrentElement, double &residual, double xiTolerance)
{
  const element_no_t nElements = this->nElementsLocal();
  VLOG(1) << "findPosition(" << point << ", elementNoLocal: " << elementNoLocal << ", ghostMeshNo: " << ghostMeshNo
    << ", startSearchInCurrentElement: " << startSearchInCurrentElement << ", xiTolerance: " << xiTolerance << ", xi: " << xi << ")";

  // set starting no to 0 if it was not given and is thus arbitrarily initialized
  if (elementNoLocal < 0 || elementNoLocal >= nElements)
    elementNoLocal = 0;

  // variables to store the best found element so far
  element_no_t elementNoBest = 0;
  std::array<double,MeshType::dim()> xiBest;
  double excessivityScoreBest = std::numeric_limits<double>::max();
  double residualBest = 0;
  bool elementFound = false;

  if (startSearchInCurrentElement)
  {
    VLOG(1) << "findPosition: startSearchInCurrentElement";

    FunctionSpaceStructuredFindPositionBase<MeshType,BasisFunctionType> *functionSpace = this;

    if (ghostMeshNo != -1)
    {
      functionSpace = ghostMesh_[ghostMeshNo].get();

      // if no such ghost mesh is stored, fall back to own mesh
      if (functionSpace == nullptr)
      {
        ghostMeshNo = -1;
        functionSpace = this;
      }
    }

    assert(functionSpace);

    // check if point is already in current element
    if (functionSpace->pointIsInElement(point, elementNoLocal, xi, residual, xiTolerance))
    {
      elementFound = true;

      // debugging output
      if (VLOG_IS_ON(2) && MeshType::dim() == 3)
      {
        // check for element size
        const int nDofsPerElement = FunctionSpace<MeshType,BasisFunctionType>::nDofsPerElement();
        std::array<Vec3,nDofsPerElement> elementalNodePositions;
        this->geometryField().getElementValues(elementNoLocal, elementalNodePositions);

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

        VLOG(1) << "point " << point << " is in element " << elementNoLocal << ", which has "
          << "bounding box x: [" << xMin << "," << xMax << "], y: [" << yMin << "," << yMax << "], z: [" << zMin << "," << zMax << "]";
      }
      
      // check if point is really inside the element by a tigher tolerance
      double excessivityScore = 0;      // lower means more inside the element, <= 0 equals really totally inside the element, then (0 < xi < 1)
      for (int i = 0; i < MeshType::dim(); i++)
      {
        excessivityScore = std::max({excessivityScore, xi[i] - 1.0, 0.0 - xi[i]});
      }

      
      // if the point is really inside the element even with the tight tolerance, return true,
      // otherwise look in neighbouring elements for a better fit
      if (excessivityScore < 1e-12)
      {
        VLOG(1) << "findPosition: pointIsInElement returned true, found at xi=" << xi << ", elementNo: " << elementNoLocal << ", excessivityScore=" << excessivityScore << ", use it";
        return true;
      }
      else 
      {
        VLOG(1) << "findPosition: pointIsInElement returned true, found at xi=" << xi << ", elementNo: " << elementNoLocal << ", excessivityScore=" << excessivityScore << ", save and check neighbouring elements";
        
        if (excessivityScore < excessivityScoreBest)
        {
          // save element as the best one so far, but also check neighbouring elements
          elementNoBest = elementNoLocal;
          xiBest = xi;
          residualBest = residual;
          excessivityScoreBest = excessivityScore;
          VLOG(1) << "findPosition: stored element " << elementNoBest << ", xi " << xiBest << ", residual: " << residualBest << ", score: " << excessivityScoreBest;
        }
      }
    }

    // point is not in current element, consider the neighbouring elements and ghost meshes

    VLOG(1) << "point is not in current element " << elementNoLocal << ", now check neighbouring elements";

    // set the neighbouring element nos, also considering ghost meshes
    if (this->checkNeighbouringElements(point, elementNoLocal, ghostMeshNo, xi, residual))
    {
      elementFound = true;

      // determine excessivity score of point found in neighbouring elements 
      double excessivityScore = 0;      // lower means more inside the element, <= 0 equals really totally inside the element, then (0 < xi < 1)
      for (int i = 0; i < MeshType::dim(); i++)
      {
        excessivityScore = std::max({excessivityScore, xi[i] - 1.0, 0.0 - xi[i]});
      }

      VLOG(1) << "findPosition: checkNeighbouringElements returned true, found at xi=" << xi << ", elementNo: " << elementNoLocal << ", excessivityScore=" << excessivityScore;

      if (excessivityScore < excessivityScoreBest)
      {
        VLOG(1) << "findPosition: use the one found (of neighbour element)";

        // the element of the neighbour was better, use it
        return true;
      }
      else 
      {
        VLOG(1) << "findPosition: use previously found, that is better";
        
        // use previously found element which is better
        xi = xiBest;
        residual = residualBest;
        elementNoLocal = elementNoBest;
        return true;
      }
    }
    else
    {
      VLOG(1) << "point was also not found among neighbouring elements (including ghost meshes), xi: " << xi << ", now check all elements";
    }
  }

  if (elementFound)
  {
    elementNoLocal = elementNoBest;
    xi = xiBest;
    residual = residualBest;
    
    VLOG(1) << "findPosition: element was found earlier with xi=" << xi << ", elementNo: " << elementNoLocal << ", excessivityScore=" << excessivityScoreBest << ", use it.";

    return true;
  }

  // search among all elements

  // look in every element, starting at elementNoLocal-2
  element_no_t elementNoLocalStart = (elementNoLocal - 2 + nElements) % nElements;
  element_no_t elementNoLocalEnd = elementNoLocalStart;

  VLOG(1) << "elementNoLocalStart: " << elementNoLocalStart << ", elementNoLocalEnd: " << elementNoLocalEnd << ", nElements: " << nElements;
  if (this->dim() == 3)
    VLOG(1) << "(" << this->meshPartition_->nElementsLocal(0) << "x" << this->meshPartition_->nElementsLocal(1) << "x" << this->meshPartition_->nElementsLocal(2) << ")";
  if (this->dim() == 2)
    VLOG(1) << "(" << this->meshPartition_->nElementsLocal(0) << "x" << this->meshPartition_->nElementsLocal(1) << ")";

  int nElementsVisited = 0;
  for (element_no_t currentElementNo = elementNoLocalStart; currentElementNo != elementNoLocalEnd || nElementsVisited == 0; currentElementNo++, nElementsVisited++)
  {
    nElementsVisited++;

    // restart with element 0
    if (currentElementNo == nElements)
    {
      currentElementNo = 0;
      if (elementNoLocalEnd == currentElementNo)
        break;
    }

    VLOG(1) << "check element " << currentElementNo;

    // check if point is already in current element
    if (this->pointIsInElement(point, currentElementNo, xi, residual, xiTolerance))
    {
      elementFound = true;

      if (startSearchInCurrentElement)
      {
#ifndef NDEBUG
        LOG(WARNING) << "Mesh \"" << this->meshName() << "\": Could not find element that contains point " << point
          << " in neighbourhood of element " << elementNoLocal
          << (ghostMeshNo != -1? std::string(" in ghost mesh ")+Mesh::getString((Mesh::face_t)ghostMeshNo) : "")
          << ", tested all elements (no ghost elements) and found element " << currentElementNo << ". "
          << "This can happen if the elements lies on the border of the higher dimensional element, e.g. if a fiber lies on the outer border of the 3D muscle mesh.";
#endif
      }

      // check if point is really inside the element by a tigher tolerance
      double excessivityScore = 0;      // lower means more inside the element, <= 0 equals really totally inside the element, then (0 < xi < 1)
      for (int i = 0; i < MeshType::dim(); i++)
      {
        excessivityScore = std::max({excessivityScore, xi[i] - 1.0, 0.0 - xi[i]});
      }
      
      // if the point is really inside the element even with the tight tolerance, return true,
      // otherwise look in neighbouring elements for a better fit
      if (excessivityScore < 1e-12)
      {
        VLOG(1) << "findPosition, checking all elements: pointIsInElement returned true, found at xi=" << xi << ", elementNo: " << elementNoLocal << ", excessivityScore=" << excessivityScore << ", use it";
        elementNoLocal = currentElementNo;
        return true;
      }
      else 
      {
        VLOG(1) << "findPosition, checking all elements: pointIsInElement returned true, found at xi=" << xi << ", elementNo: " << currentElementNo << ", excessivityScore=" << excessivityScore << ", save";
        // save element as the best one so far, but also check neighbouring elements
        if (excessivityScore < excessivityScoreBest)
        {
          elementNoBest = currentElementNo;
          xiBest = xi;
          residualBest = residual;
          excessivityScoreBest = excessivityScore;
          ghostMeshNo = -1;   // not a ghost mesh
          VLOG(1) << "findPosition, checking all elements: stored element " << elementNoBest << ", xi " << xiBest << ", residual: " << residualBest << ", score: " << excessivityScoreBest;
        }
      }
    }
  }

  if (elementFound)
  {
    elementNoLocal = elementNoBest;
    xi = xiBest;
    residual = residualBest;
    
    VLOG(1) << "findPosition: element was found earlier with xi=" << xi << ", elementNo: " << elementNoLocal << ", excessivityScore=" << excessivityScoreBest << ", use it.";

    return true;
  }
  
  // if point was still not found, search in ghost meshes
  for (int face = (int)Mesh::face_t::face0Minus; face <= (int)Mesh::face_t::face2Plus; face++)
  {
    VLOG(3) << "consider ghost mesh " << Mesh::getString((Mesh::face_t)face);
    if (ghostMesh_[face] != nullptr)
    {
      VLOG(3) << "   ghost mesh " << Mesh::getString((Mesh::face_t)face) << " is set";
      if (ghostMesh_[face]->findPosition(point, elementNoLocal, ghostMeshNo, xi, false, residual))
      {
        VLOG(3) << "   point found in ghost mesh " << Mesh::getString((Mesh::face_t)face) << ", element " << elementNoLocal << ", xi " << xi;
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

} // namespace
