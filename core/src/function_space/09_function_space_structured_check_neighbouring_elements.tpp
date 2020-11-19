#include "function_space/09_function_space_structured_check_neighbouring_elements.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"
#include "mesh/face_t.h"

namespace FunctionSpace
{

// 1D
template<typename MeshType, typename BasisFunctionType>
bool FunctionSpaceStructuredCheckNeighbouringElements<MeshType,BasisFunctionType,Mesh::isDim<1,MeshType>>::
checkNeighbouringElements(const Vec3 &point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,MeshType::dim()> &xi, double &residual, double xiTolerance)
{
  VLOG(1) << "checkNeighbouringElements<1D>(elementNo = " << elementNo << ", ghostMeshNo = " << ghostMeshNo << ", initial xi = " << xi;

  // variables to store the best found element so far
  element_no_t elementNoBest = 0;
  std::array<double,MeshType::dim()> xiBest;
  double excessivityScoreBest = std::numeric_limits<double>::max();
  double residualBest;
  bool elementFound = false;

  ghostMeshNo = -1;
  for (element_no_t neighbourElementNo = elementNo - 1; neighbourElementNo != elementNo + 2; neighbourElementNo++)
  {
    if (0 <= neighbourElementNo && neighbourElementNo < this->meshPartition_->nElementsLocal(0))
    {
      if (this->pointIsInElement(point, neighbourElementNo, xi, residual, xiTolerance))
      {
        elementFound = true;

        // check if found point is really inside the element by a tighter tolerance
        // lower means more inside the element, <= 0 equals really totally inside the element, then (0 < xi < 1)
        double excessivityScore = std::max(xi[0] - 1.0, 0.0 - xi[0]);

        // if the point is really inside the element even with the tight tolerance, return true,
        // otherwise look in neighbouring elements for a better fit
        if (excessivityScore < 1e-12)
        {
          elementNo = neighbourElementNo;
          return true;
        }
        else 
        {
          // if the found element is the best so far
          if (excessivityScore < excessivityScoreBest)
          {
            // save element as the best one so far, but also check neighbouring elements
            elementNoBest = neighbourElementNo;
            xiBest = xi;
            residualBest = residual;
            excessivityScoreBest = excessivityScore;
          }
        }
      }
    }
  }

  // if at least one element was found that contains the point, use it and return true
  if (elementFound)
  {
    elementNo = elementNoBest;
    xi = xiBest;
    residual = residualBest;
    return true;
  }

  return false;
}

// 2D
template<typename MeshType, typename BasisFunctionType>
bool FunctionSpaceStructuredCheckNeighbouringElements<MeshType,BasisFunctionType,Mesh::isDim<2,MeshType>>::
checkNeighbouringElements(const Vec3 &point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,MeshType::dim()> &xi, double &residual, double xiTolerance)
{
  VLOG(1) << "checkNeighbouringElements<2D>(elementNo = " << elementNo << ", ghostMeshNo = " << ghostMeshNo << ", initial xi = " << xi;

  const int D = 2;

  // define the order in which the neighbors are considered
  static std::array<int,3> xOffset;
  static std::array<int,3> yOffset;

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

  // get local coordinates of element
  std::array<int,2> coordinatesLocal;
  if (ghostMeshNo == -1)
  {
    coordinatesLocal = this->meshPartition_->getElementCoordinatesLocal(elementNo);
  }
  else if (ghostMeshNo == (int)Mesh::face_t::face0Minus)
  {
    assert(this->ghostMesh_[ghostMeshNo]);
    coordinatesLocal = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementCoordinatesLocal(elementNo);
    coordinatesLocal[0] -= this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0);
  }
  else if (ghostMeshNo == (int)Mesh::face_t::face0Plus)
  {
    assert(this->ghostMesh_[ghostMeshNo]);
    coordinatesLocal = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementCoordinatesLocal(elementNo);
    coordinatesLocal[0] += this->meshPartition_->nElementsLocal(0);
  }
  else if (ghostMeshNo == (int)Mesh::face_t::face1Minus)
  {
    assert(this->ghostMesh_[ghostMeshNo]);
    coordinatesLocal = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementCoordinatesLocal(elementNo);
    coordinatesLocal[1] -= this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1);
  }
  else if (ghostMeshNo == (int)Mesh::face_t::face1Plus)
  {
    assert(this->ghostMesh_[ghostMeshNo]);
    coordinatesLocal = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementCoordinatesLocal(elementNo);
    coordinatesLocal[1] += this->meshPartition_->nElementsLocal(1);
  }
  else if (ghostMeshNo == (int)Mesh::face_or_edge_t::edge0Minus1Minus)  // bottom left
  {
    assert(this->ghostMesh_[ghostMeshNo]);
    coordinatesLocal = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementCoordinatesLocal(elementNo);
    coordinatesLocal[0] -= this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0);
    coordinatesLocal[1] -= this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1);
  }
  else if (ghostMeshNo == (int)Mesh::face_or_edge_t::edge0Plus1Minus)  // bottom right
  {
    assert(this->ghostMesh_[ghostMeshNo]);
    coordinatesLocal = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementCoordinatesLocal(elementNo);
    coordinatesLocal[0] += this->meshPartition_->nElementsLocal(0);
    coordinatesLocal[1] -= this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1);
  }
  else if (ghostMeshNo == (int)Mesh::face_or_edge_t::edge0Minus1Plus)  // top left
  {
    assert(this->ghostMesh_[ghostMeshNo]);
    coordinatesLocal = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementCoordinatesLocal(elementNo);
    coordinatesLocal[0] -= this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0);
    coordinatesLocal[1] += this->meshPartition_->nElementsLocal(1);
  }
  else if (ghostMeshNo == (int)Mesh::face_or_edge_t::edge0Plus1Plus)  // top right
  {
    assert(this->ghostMesh_[ghostMeshNo]);
    coordinatesLocal = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementCoordinatesLocal(elementNo);
    coordinatesLocal[0] += this->meshPartition_->nElementsLocal(0);
    coordinatesLocal[1] += this->meshPartition_->nElementsLocal(1);
  }

  VLOG(1) << "nElementsLocal: [" << this->meshPartition_->nElementsLocal(0) << "," << this->meshPartition_->nElementsLocal(1)
    << "] coordinatesLocal: " << coordinatesLocal
    << " interation y in " << yOffset << " x in " << xOffset << "";

  //debugOutputGhostMeshSet();

  // variables to store the best found element so far
  element_no_t elementNoBest = 0;
  int ghostMeshNoBest = 0;
  std::array<double,MeshType::dim()> xiBest;
  double excessivityScoreBest = std::numeric_limits<double>::max();
  bool elementFound = false;
  double residualBest;

  std::array<int,2> neighbourCoordinatesLocal;
  element_no_t neighbourElementNo;

  FunctionSpaceStructuredFindPositionBase<MeshType,BasisFunctionType> *functionSpace = this;

  for (int yIndex = 0; yIndex != 3; yIndex++)
  {
    int y = yOffset[yIndex];
    for (int xIndex = 0; xIndex != 3; xIndex++)
    {
      int x = xOffset[xIndex];

      // do not handle the center point
      if (x == 0 && y == 0)
      {
        continue;
      }

      neighbourElementNo = -1;   // set to invalid

      // determine coordinates of neighbor element that is currently examined
      neighbourCoordinatesLocal = coordinatesLocal + std::array<int,2>({x,y});

      VLOG(1) << "(x,y) = (" << x << "," << y << "), neighbourCoordinatesLocal: " << neighbourCoordinatesLocal;

      functionSpace = this;

      if (neighbourCoordinatesLocal[0] < 0
          && neighbourCoordinatesLocal[1] < 0)  // if at bottom left corner
      {
        if (this->ghostMesh_[(int)Mesh::face_or_edge_t::edge0Minus1Minus] != nullptr)
        {
          ghostMeshNo = (int)Mesh::face_or_edge_t::edge0Minus1Minus;
          functionSpace = this->ghostMesh_[ghostMeshNo].get();
          std::array<int,2> ghostMeshCoordinates(neighbourCoordinatesLocal);
          ghostMeshCoordinates[0] += this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0);
          ghostMeshCoordinates[1] += this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1);
          neighbourElementNo = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);

          VLOG(1) << "0-1- neighbourElementNo: " << neighbourElementNo << ", ghostMeshCoordinates: " << ghostMeshCoordinates
            << " / (" << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(2) << ")";
        }
        else
        {
          // at this location outside of the normal subdomain no ghost mesh was specified
          continue;
        }
      }
      else if (neighbourCoordinatesLocal[0] >= this->meshPartition_->nElementsLocal(0)
               && neighbourCoordinatesLocal[1] < 0)  // if at bottom right corner
      {
        if (this->ghostMesh_[(int)Mesh::face_or_edge_t::edge0Plus1Minus] != nullptr)
        {
          ghostMeshNo = (int)Mesh::face_or_edge_t::edge0Plus1Minus;
          functionSpace = this->ghostMesh_[ghostMeshNo].get();
          std::array<int,2> ghostMeshCoordinates(neighbourCoordinatesLocal);
          ghostMeshCoordinates[0] -= this->meshPartition_->nElementsLocal(0);
          ghostMeshCoordinates[1] += this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1);
          neighbourElementNo = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);

          VLOG(1) << "0+1- neighbourElementNo: " << neighbourElementNo << ", ghostMeshCoordinates: " << ghostMeshCoordinates
            << " / (" << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(2) << ")";
        }
        else
        {
          // at this location outside of the normal subdomain no ghost mesh was specified
          continue;
        }
      }
      else if (neighbourCoordinatesLocal[0] < 0
          && neighbourCoordinatesLocal[1] >= this->meshPartition_->nElementsLocal(1))  // if at top left corner
      {
        if (this->ghostMesh_[(int)Mesh::face_or_edge_t::edge0Minus1Plus] != nullptr)
        {
          ghostMeshNo = (int)Mesh::face_or_edge_t::edge0Minus1Plus;
          functionSpace = this->ghostMesh_[ghostMeshNo].get();
          std::array<int,2> ghostMeshCoordinates(neighbourCoordinatesLocal);
          ghostMeshCoordinates[0] += this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0);
          ghostMeshCoordinates[1] -= this->meshPartition_->nElementsLocal(1);
          neighbourElementNo = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);

          VLOG(1) << "0-1+ neighbourElementNo: " << neighbourElementNo << ", ghostMeshCoordinates: " << ghostMeshCoordinates
            << " / (" << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(2) << ")";
        }
        else
        {
          // at this location outside of the normal subdomain no ghost mesh was specified
          continue;
        }
      }
      else if (neighbourCoordinatesLocal[0] >= this->meshPartition_->nElementsLocal(0)
          && neighbourCoordinatesLocal[1] >= this->meshPartition_->nElementsLocal(1))  // if at top right corner
      {
        if (this->ghostMesh_[(int)Mesh::face_or_edge_t::edge0Plus1Plus] != nullptr)
        {
          ghostMeshNo = (int)Mesh::face_or_edge_t::edge0Plus1Plus;
          functionSpace = this->ghostMesh_[ghostMeshNo].get();
          std::array<int,2> ghostMeshCoordinates(neighbourCoordinatesLocal);
          ghostMeshCoordinates[0] -= this->meshPartition_->nElementsLocal(0);
          ghostMeshCoordinates[1] -= this->meshPartition_->nElementsLocal(1);
          neighbourElementNo = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);

          VLOG(1) << "0+1+ neighbourElementNo: " << neighbourElementNo << ", ghostMeshCoordinates: " << ghostMeshCoordinates
            << " / (" << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(2) << ")";
        }
        else
        {
          // at this location outside of the normal subdomain no ghost mesh was specified
          continue;
        }
      }
      else if (neighbourCoordinatesLocal[0] < 0)  // if at left boundary
      {
        if (this->ghostMesh_[(int)Mesh::face_t::face0Minus] != nullptr)
        {
          ghostMeshNo = (int)Mesh::face_t::face0Minus;
          functionSpace = this->ghostMesh_[ghostMeshNo].get();
          std::array<int,2> ghostMeshCoordinates(neighbourCoordinatesLocal);
          ghostMeshCoordinates[0] += this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0);
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
      else if (neighbourCoordinatesLocal[0] >= this->meshPartition_->nElementsLocal(0))  // if at right boundary
      {
        if (this->ghostMesh_[(int)Mesh::face_t::face0Plus] != nullptr)
        {
          ghostMeshNo = (int)Mesh::face_t::face0Plus;
          functionSpace = this->ghostMesh_[ghostMeshNo].get();
          std::array<int,2> ghostMeshCoordinates(neighbourCoordinatesLocal);
          ghostMeshCoordinates[0] -= this->meshPartition_->nElementsLocal(0);
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
      else if (neighbourCoordinatesLocal[1] < 0)  // if at bottom boundary
      {
        if (this->ghostMesh_[(int)Mesh::face_t::face1Minus] != nullptr)
        {
          ghostMeshNo = (int)Mesh::face_t::face1Minus;
          functionSpace = this->ghostMesh_[ghostMeshNo].get();
          std::array<int,2> ghostMeshCoordinates(neighbourCoordinatesLocal);
          ghostMeshCoordinates[1] += this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1);
          neighbourElementNo = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);

          VLOG(1) << "1- neighbourElementNo: " << neighbourElementNo << ", ghostMeshCoordinates: " << ghostMeshCoordinates
            << " / (" << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(2) << ")";
        }
      }
      else if (neighbourCoordinatesLocal[1] >= this->meshPartition_->nElementsLocal(1))  // if at top boundary
      {
        if (this->ghostMesh_[(int)Mesh::face_t::face1Plus] != nullptr)
        {
          ghostMeshNo = (int)Mesh::face_t::face1Plus;
          functionSpace = this->ghostMesh_[ghostMeshNo].get();
          std::array<int,2> ghostMeshCoordinates(neighbourCoordinatesLocal);
          ghostMeshCoordinates[1] -= this->meshPartition_->nElementsLocal(1);
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
      else if (0 <= neighbourCoordinatesLocal[0] && neighbourCoordinatesLocal[0] < this->meshPartition_->nElementsLocal(0)
        && 0 <= neighbourCoordinatesLocal[1] && neighbourCoordinatesLocal[1] < this->meshPartition_->nElementsLocal(1))
      {
        neighbourElementNo = this->meshPartition_->getElementNoLocal(neighbourCoordinatesLocal);
        functionSpace = this;
        ghostMeshNo = -1;
        VLOG(1) << "normal, neighbourElementNo: " << neighbourElementNo << " nElementsLocal: " << this->nElementsLocal() << "=" << functionSpace->nElementsLocal();
      }

      // check if point is in current neighbour element
      if (neighbourElementNo != -1)
      {
        if (functionSpace->pointIsInElement(point, neighbourElementNo, xi, residual, xiTolerance))
        {
          elementFound = true;

          // check if found point is really inside the element by a tighter tolerance
          double excessivityScore = 0;      // lower means more inside the element, <= 0 equals really totally inside the element, then (0 < xi < 1)
          for (int i = 0; i < D; i++)
          {
            excessivityScore = std::max({excessivityScore, xi[i] - 1.0, 0.0 - xi[i]});
          }

          // if the point is really inside the element even with the tight tolerance, return true,
          // otherwise continue looking in elements for a better fit
          if (excessivityScore < 1e-12)
          {
            elementNo = neighbourElementNo;
            return true;
          }
          else 
          {
            // if the found element is the best so far
            if (excessivityScore < excessivityScoreBest)
            {
              // save element as the best one so far, but also check neighbouring elements
              elementNoBest = neighbourElementNo;
              xiBest = xi;
              ghostMeshNoBest = ghostMeshNo;
              residualBest = residual;
              excessivityScoreBest = excessivityScore;
            }
          }
        }
        else
        {
          VLOG(1) << "  point is not in element " << neighbourElementNo << " (xi: " << xi << ").";
        }
      }
    }  // x
  }   // y

  // if at least one element was found that contains the point, use it and return true
  if (elementFound)
  {
    elementNo = elementNoBest;
    xi = xiBest;
    ghostMeshNo = ghostMeshNoBest;
    residual = residualBest;
    return true;
  }

  return false;
}

// 3D
template<typename MeshType, typename BasisFunctionType>
bool FunctionSpaceStructuredCheckNeighbouringElements<MeshType,BasisFunctionType,Mesh::isDim<3,MeshType>>::
checkNeighbouringElements(const Vec3 &point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,MeshType::dim()> &xi, double &residual, double xiTolerance)
{
  LOG(DEBUG) << "  " << "checkNeighbouringElements(elementNo = " << elementNo << ", ghostMeshNo = " << ghostMeshNo << ", initial xi = " << xi;

  const int D = 3;

  // define the order in which the neighbors are considered
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
  else if (ghostMeshNo == (int)Mesh::face_t::face0Minus)  // left
  {
    assert(this->ghostMesh_[ghostMeshNo]);
    coordinatesLocal = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementCoordinatesLocal(elementNo);
    coordinatesLocal[0] -= -this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0);
  }
  else if (ghostMeshNo == (int)Mesh::face_t::face0Plus)   // right
  {
    assert(this->ghostMesh_[ghostMeshNo]);
    coordinatesLocal = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementCoordinatesLocal(elementNo);
    coordinatesLocal[0] += this->meshPartition_->nElementsLocal(0);
  }
  else if (ghostMeshNo == (int)Mesh::face_t::face1Minus)  // front
  {
    assert(this->ghostMesh_[ghostMeshNo]);
    coordinatesLocal = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementCoordinatesLocal(elementNo);
    coordinatesLocal[1] -= this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1);
  }
  else if (ghostMeshNo == (int)Mesh::face_t::face1Plus)   // back
  {
    assert(this->ghostMesh_[ghostMeshNo]);
    coordinatesLocal = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementCoordinatesLocal(elementNo);
    coordinatesLocal[1] += this->meshPartition_->nElementsLocal(1);
  }
  else if (ghostMeshNo == (int)Mesh::face_t::face2Minus)  // bottom
  {
    assert(this->ghostMesh_[ghostMeshNo]);
    coordinatesLocal = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementCoordinatesLocal(elementNo);
    coordinatesLocal[2] -= this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(2);
  }
  else if (ghostMeshNo == (int)Mesh::face_t::face2Plus)   // top
  {
    assert(this->ghostMesh_[ghostMeshNo]);
    coordinatesLocal = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementCoordinatesLocal(elementNo);
    coordinatesLocal[2] += this->meshPartition_->nElementsLocal(2);
  }
  else if (ghostMeshNo == (int)Mesh::face_or_edge_t::edge0Minus1Minus)  // front left
  {
    assert(this->ghostMesh_[ghostMeshNo]);
    coordinatesLocal = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementCoordinatesLocal(elementNo);
    coordinatesLocal[0] -= this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0);
    coordinatesLocal[1] -= this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1);
  }
  else if (ghostMeshNo == (int)Mesh::face_or_edge_t::edge0Plus1Minus)  // front right
  {
    assert(this->ghostMesh_[ghostMeshNo]);
    coordinatesLocal = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementCoordinatesLocal(elementNo);
    coordinatesLocal[0] += this->meshPartition_->nElementsLocal(0);
    coordinatesLocal[1] -= this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1);
  }
  else if (ghostMeshNo == (int)Mesh::face_or_edge_t::edge0Minus1Plus)  // back left
  {
    assert(this->ghostMesh_[ghostMeshNo]);
    coordinatesLocal = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementCoordinatesLocal(elementNo);
    coordinatesLocal[0] -= this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0);
    coordinatesLocal[1] += this->meshPartition_->nElementsLocal(1);
  }
  else if (ghostMeshNo == (int)Mesh::face_or_edge_t::edge0Plus1Plus)  // back right
  {
    assert(this->ghostMesh_[ghostMeshNo]);
    coordinatesLocal = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementCoordinatesLocal(elementNo);
    coordinatesLocal[0] += this->meshPartition_->nElementsLocal(0);
    coordinatesLocal[1] += this->meshPartition_->nElementsLocal(1);
  }

  VLOG(2) << "nElementsLocal: [" << this->meshPartition_->nElementsLocal(0) << "," << this->meshPartition_->nElementsLocal(1)
    << "," << this->meshPartition_->nElementsLocal(2) << "] coordinatesLocal: " << coordinatesLocal
    << " interation z in " << zOffset << ", y in " << yOffset << " x in " << xOffset << "";

  //debugOutputGhostMeshSet();

  // variables to store the best found element so far
  element_no_t elementNoBest = 0;
  int ghostMeshNoBest = 0;
  std::array<double,MeshType::dim()> xiBest; 
  double residualBest;
  double excessivityScoreBest = std::numeric_limits<double>::max();
  LOG(DEBUG) << "intiialize excessivityScoreBest: " << excessivityScoreBest;

  bool elementFound = false;

  std::array<int,3> neighbourCoordinatesLocal;
  element_no_t neighbourElementNo;

  FunctionSpaceStructuredFindPositionBase<MeshType,BasisFunctionType> *functionSpace = this;

  // loop over neighbors of current element
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

        // determine coordinates of neighbor element that is currently examined
        neighbourCoordinatesLocal = coordinatesLocal + std::array<int,3>({x,y,z});

        VLOG(2) << "(x,y,z) = (" << x << "," << y << "," << z << "), neighbourCoordinatesLocal: " << neighbourCoordinatesLocal;
        functionSpace = this;

        if (neighbourCoordinatesLocal[0] < 0
            && neighbourCoordinatesLocal[1] < 0)  // if at front left corner
        {
          if (this->ghostMesh_[(int)Mesh::face_or_edge_t::edge0Minus1Minus] != nullptr)
          {
            ghostMeshNo = (int)Mesh::face_or_edge_t::edge0Minus1Minus;
            functionSpace = this->ghostMesh_[ghostMeshNo].get();
            std::array<int,3> ghostMeshCoordinates(neighbourCoordinatesLocal);
            ghostMeshCoordinates[0] += this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0);
            ghostMeshCoordinates[1] += this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1);
            neighbourElementNo = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);

            VLOG(1) << "0-1- neighbourElementNo: " << neighbourElementNo << ", ghostMeshCoordinates: " << ghostMeshCoordinates
              << " / (" << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(2) << ")";
          }
          else
          {
            // at this location outside of the normal subdomain no ghost mesh was specified
            continue;
          }
        }
        else if (neighbourCoordinatesLocal[0] >= this->meshPartition_->nElementsLocal(0)
                && neighbourCoordinatesLocal[1] < 0)  // if at front right corner
        {
          if (this->ghostMesh_[(int)Mesh::face_or_edge_t::edge0Plus1Minus] != nullptr)
          {
            ghostMeshNo = (int)Mesh::face_or_edge_t::edge0Plus1Minus;
            functionSpace = this->ghostMesh_[ghostMeshNo].get();
            std::array<int,3> ghostMeshCoordinates(neighbourCoordinatesLocal);
            ghostMeshCoordinates[0] -= this->meshPartition_->nElementsLocal(0);
            ghostMeshCoordinates[1] += this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1);
            neighbourElementNo = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);

            VLOG(1) << "0+1- neighbourElementNo: " << neighbourElementNo << ", ghostMeshCoordinates: " << ghostMeshCoordinates
              << " / (" << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(2) << ")";
          }
          else
          {
            // at this location outside of the normal subdomain no ghost mesh was specified
            continue;
          }
        }
        else if (neighbourCoordinatesLocal[0] < 0
            && neighbourCoordinatesLocal[1] >= this->meshPartition_->nElementsLocal(1))  // if at back left corner
        {
          if (this->ghostMesh_[(int)Mesh::face_or_edge_t::edge0Minus1Plus] != nullptr)
          {
            ghostMeshNo = (int)Mesh::face_or_edge_t::edge0Minus1Plus;
            functionSpace = this->ghostMesh_[ghostMeshNo].get();
            std::array<int,3> ghostMeshCoordinates(neighbourCoordinatesLocal);
            ghostMeshCoordinates[0] += this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0);
            ghostMeshCoordinates[1] -= this->meshPartition_->nElementsLocal(1);
            neighbourElementNo = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);

            VLOG(1) << "0-1+ neighbourElementNo: " << neighbourElementNo << ", ghostMeshCoordinates: " << ghostMeshCoordinates
              << " / (" << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(2) << ")";
          }
          else
          {
            // at this location outside of the normal subdomain no ghost mesh was specified
            continue;
          }
        }
        else if (neighbourCoordinatesLocal[0] >= this->meshPartition_->nElementsLocal(0)
            && neighbourCoordinatesLocal[1] >= this->meshPartition_->nElementsLocal(1))  // if at back right corner
        {
          if (this->ghostMesh_[(int)Mesh::face_or_edge_t::edge0Plus1Plus] != nullptr)
          {
            ghostMeshNo = (int)Mesh::face_or_edge_t::edge0Plus1Plus;
            functionSpace = this->ghostMesh_[ghostMeshNo].get();
            std::array<int,3> ghostMeshCoordinates(neighbourCoordinatesLocal);
            ghostMeshCoordinates[0] -= this->meshPartition_->nElementsLocal(0);
            ghostMeshCoordinates[1] -= this->meshPartition_->nElementsLocal(1);
            neighbourElementNo = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);

            VLOG(1) << "0+1+ neighbourElementNo: " << neighbourElementNo << ", ghostMeshCoordinates: " << ghostMeshCoordinates
              << " / (" << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(2) << ")";
          }
          else
          {
            // at this location outside of the normal subdomain no ghost mesh was specified
            continue;
          }
        }
        else if (neighbourCoordinatesLocal[0] < 0)  // if at left boundary
        {
          if (this->ghostMesh_[(int)Mesh::face_t::face0Minus] != nullptr)
          {
            ghostMeshNo = (int)Mesh::face_t::face0Minus;
            functionSpace = this->ghostMesh_[ghostMeshNo].get();
            std::array<int,3> ghostMeshCoordinates(neighbourCoordinatesLocal);
            ghostMeshCoordinates[0] += this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0);
            neighbourElementNo = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);

            VLOG(2) << "0- neighbourElementNo: " << neighbourElementNo << ", ghostMeshCoordinates: " << ghostMeshCoordinates
              << " / (" << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(2) << ")";
          }
          else
          {
            // at this position outside of the normal subdomain no ghost mesh was specified
            continue;
          }
        }
        else if (neighbourCoordinatesLocal[0] >= this->meshPartition_->nElementsLocal(0))  // if at right boundary
        {
          if (this->ghostMesh_[(int)Mesh::face_t::face0Plus] != nullptr)
          {
            ghostMeshNo = (int)Mesh::face_t::face0Plus;
            functionSpace = this->ghostMesh_[ghostMeshNo].get();
            std::array<int,3> ghostMeshCoordinates(neighbourCoordinatesLocal);
            ghostMeshCoordinates[0] -= this->meshPartition_->nElementsLocal(0);
            neighbourElementNo = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);

            VLOG(2) << "0+ neighbourElementNo: " << neighbourElementNo << ", ghostMeshCoordinates: " << ghostMeshCoordinates
              << " / (" << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(2) << ")";
          }
          else
          {
            // at this position outside of the normal subdomain no ghost mesh was specified
            continue;
          }
        }
        else if (neighbourCoordinatesLocal[1] < 0)  // if at front boundary
        {
          if (this->ghostMesh_[(int)Mesh::face_t::face1Minus] != nullptr)
          {
            ghostMeshNo = (int)Mesh::face_t::face1Minus;
            functionSpace = this->ghostMesh_[ghostMeshNo].get();
            std::array<int,3> ghostMeshCoordinates(neighbourCoordinatesLocal);
            ghostMeshCoordinates[1] += this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1);
            neighbourElementNo = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);

            VLOG(2) << "1- neighbourElementNo: " << neighbourElementNo << ", ghostMeshCoordinates: " << ghostMeshCoordinates
              << " / (" << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(2) << ")";
          }
        }
        else if (neighbourCoordinatesLocal[1] >= this->meshPartition_->nElementsLocal(1))  // if at back boundary
        {
          if (this->ghostMesh_[(int)Mesh::face_t::face1Plus] != nullptr)
          {
            ghostMeshNo = (int)Mesh::face_t::face1Plus;
            functionSpace = this->ghostMesh_[ghostMeshNo].get();
            std::array<int,3> ghostMeshCoordinates(neighbourCoordinatesLocal);
            ghostMeshCoordinates[1] -= this->meshPartition_->nElementsLocal(1);
            neighbourElementNo = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);

            VLOG(2) << "1+ neighbourElementNo: " << neighbourElementNo << ", ghostMeshCoordinates: " << ghostMeshCoordinates
              << " / (" << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(2) << ")";
          }
          else
          {
            // at this position outside of the normal subdomain no ghost mesh was specified
            continue;
          }
        }
        else if (neighbourCoordinatesLocal[2] < 0)  // if at bottom boundary
        {
          if (this->ghostMesh_[(int)Mesh::face_t::face2Minus] != nullptr)
          {
            ghostMeshNo = (int)Mesh::face_t::face2Minus;
            functionSpace = this->ghostMesh_[ghostMeshNo].get();
            std::array<int,3> ghostMeshCoordinates(neighbourCoordinatesLocal);
            ghostMeshCoordinates[2] += this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(2);
            neighbourElementNo = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);

            VLOG(2) << "2- neighbourElementNo: " << neighbourElementNo << ", ghostMeshCoordinates: " << ghostMeshCoordinates
              << " / (" << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(2) << ")";
          }
          else
          {
            // at this position outside of the normal subdomain no ghost mesh was specified
            continue;
          }
        }
        else if (neighbourCoordinatesLocal[2] >= this->meshPartition_->nElementsLocal(2))  // if at top boundary
        {
          if (this->ghostMesh_[(int)Mesh::face_t::face2Plus] != nullptr)
          {
            ghostMeshNo = (int)Mesh::face_t::face2Plus;
            functionSpace = this->ghostMesh_[ghostMeshNo].get();
            std::array<int,3> ghostMeshCoordinates(neighbourCoordinatesLocal);
            ghostMeshCoordinates[2] -= this->meshPartition_->nElementsLocal(2);
            neighbourElementNo = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);

            VLOG(2) << "2+ neighbourElementNo: " << neighbourElementNo << ", ghostMeshCoordinates: " << ghostMeshCoordinates
              << " / (" << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(2) << ")";
          }
          else
          {
            // at this position outside of the normal subdomain no ghost mesh was specified
            continue;
          }
        }
        else if (0 <= neighbourCoordinatesLocal[0] && neighbourCoordinatesLocal[0] < this->meshPartition_->nElementsLocal(0)
          && 0 <= neighbourCoordinatesLocal[1] && neighbourCoordinatesLocal[1] < this->meshPartition_->nElementsLocal(1)
          && 0 <= neighbourCoordinatesLocal[2] && neighbourCoordinatesLocal[2] < this->meshPartition_->nElementsLocal(2))
        {
          neighbourElementNo = this->meshPartition_->getElementNoLocal(neighbourCoordinatesLocal);
          functionSpace = this;
          ghostMeshNo = -1;
          LOG(DEBUG) << "  " << "normal, neighbourElementNo: " << neighbourElementNo << " nElementsLocal: " << this->nElementsLocal() << "=" << functionSpace->nElementsLocal();
        }

        // check if point is in current neighbour element
        if (neighbourElementNo != -1)
        {
          if (functionSpace->pointIsInElement(point, neighbourElementNo, xi, residual, xiTolerance))
          {
            elementFound = true;

            // check if found point is really inside the element by a tighter tolerance
            double excessivityScore = 0;      // lower means more inside the element, <= 0 equals really totally inside the element, then (0 < xi < 1)
            for (int i = 0; i < D; i++)
            {
              excessivityScore = std::max({excessivityScore, xi[i] - 1.0, 0.0 - xi[i]});
            }
            
            LOG(DEBUG) << "  " << "  checkNeighbouringElements: pointIsInElement returned true, found at xi=" << xi << ", elementNo: " << neighbourElementNo << ", excessivityScore=" << excessivityScore;

            // if the point is really inside the element even with the tight tolerance, return true,
            // otherwise continue looking in elements for a better fit
            if (excessivityScore < 1e-12)
            {
              elementNo = neighbourElementNo;
              return true;
            }
            else 
            {
              // if the found element is the best so far
              if (excessivityScore < excessivityScoreBest)
              {
                // save element as the best one so far, but also check neighbouring elements
                elementNoBest = neighbourElementNo;
                xiBest = xi;
                ghostMeshNoBest = ghostMeshNo;
                residualBest = residual;

                excessivityScoreBest = excessivityScore;
                LOG(DEBUG) << "set excessivityScoreBest to " << excessivityScoreBest;
              }
            }
          }
          else
          {
            LOG(DEBUG) << "  point is not in element " << neighbourElementNo << " (xi: " << xi << ", residual: " << residual << ") with xiTolerance: " << xiTolerance << ".";
          }
        }
      }  // x
    }  // y
  }   // z


  // if at least one element was found that contains the point, use it and return true
  if (elementFound)
  {

    elementNo = elementNoBest;
    xi = xiBest;
    ghostMeshNo = ghostMeshNoBest;
    residual = residualBest;
    LOG(DEBUG) << "  checkNeighbouringElements: best found was xi = " << xi << ", elementNo: " << elementNo;
    return true;
  }

  return false;
}

} // namespace
