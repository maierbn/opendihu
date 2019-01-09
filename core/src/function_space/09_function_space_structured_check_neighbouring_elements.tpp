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
checkNeighbouringElements(const Vec3 &point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,MeshType::dim()> &xi)
{
  VLOG(1) << "checkNeighbouringElements<1D>(elementNo = " << elementNo << ", ghostMeshNo = " << ghostMeshNo << ", initial xi = " << xi;

  ghostMeshNo = -1;
  for (element_no_t neighbourElementNo = elementNo - 1; neighbourElementNo != elementNo + 2; neighbourElementNo++)
  {
    if (0 <= neighbourElementNo && neighbourElementNo < this->meshPartition_->nElementsLocal(0))
    {
      if (this->pointIsInElement(point, neighbourElementNo, xi))
      {
        elementNo = neighbourElementNo;
        return true;
      }
    }
  }
  return false;
}

// 2D
template<typename MeshType, typename BasisFunctionType>
bool FunctionSpaceStructuredCheckNeighbouringElements<MeshType,BasisFunctionType,Mesh::isDim<2,MeshType>>::
checkNeighbouringElements(const Vec3 &point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,MeshType::dim()> &xi)
{
  VLOG(1) << "checkNeighbouringElements<2D>(elementNo = " << elementNo << ", ghostMeshNo = " << ghostMeshNo << ", initial xi = " << xi;

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

  VLOG(1) << "nElementsLocal: [" << this->meshPartition_->nElementsLocal(0) << "," << this->meshPartition_->nElementsLocal(1)
    << "] coordinatesLocal: " << coordinatesLocal
    << " interation y in " << yOffset << " x in " << xOffset << "";

  //debugOutputGhostMeshSet();

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

      neighbourCoordinatesLocal = coordinatesLocal + std::array<int,2>({x,y});

      VLOG(1) << "(x,y) = (" << x << "," << y << "), neighbourCoordinatesLocal: " << neighbourCoordinatesLocal;

      if (neighbourCoordinatesLocal[1] > this->meshPartition_->nElementsLocal(1) || neighbourCoordinatesLocal[1] < -1
        || neighbourCoordinatesLocal[0] > this->meshPartition_->nElementsLocal(0) || neighbourCoordinatesLocal[0] < -1)
      {
        VLOG(1) << "outside ghost layer";
        continue;
      }

      // do not consider diagonal ghost meshes, e.g. x+/y+
      int nGhostTargets = 0;
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
        if (this->ghostMesh_[(int)Mesh::face_t::face0Minus] != nullptr)
        {
          ghostMeshNo = (int)Mesh::face_t::face0Minus;
          functionSpace = this->ghostMesh_[ghostMeshNo].get();
          std::array<int,2> ghostMeshCoordinates({0,neighbourCoordinatesLocal[1]});
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
        if (this->ghostMesh_[(int)Mesh::face_t::face0Plus] != nullptr)
        {
          ghostMeshNo = (int)Mesh::face_t::face0Plus;
          functionSpace = this->ghostMesh_[ghostMeshNo].get();
          std::array<int,2> ghostMeshCoordinates({0,neighbourCoordinatesLocal[1]});
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
        if (this->ghostMesh_[(int)Mesh::face_t::face1Minus] != nullptr)
        {
          ghostMeshNo = (int)Mesh::face_t::face1Minus;
          functionSpace = this->ghostMesh_[ghostMeshNo].get();
          std::array<int,2> ghostMeshCoordinates({neighbourCoordinatesLocal[0],0});
          neighbourElementNo = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);
          VLOG(1) << "1- neighbourElementNo: " << neighbourElementNo << ", ghostMeshCoordinates: " << ghostMeshCoordinates
            << " / (" << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(2) << ")";
        }
      }
      else if (neighbourCoordinatesLocal[1] == this->meshPartition_->nElementsLocal(1))  // if at back boundary
      {
        if (this->ghostMesh_[(int)Mesh::face_t::face1Plus] != nullptr)
        {
          ghostMeshNo = (int)Mesh::face_t::face1Plus;
          functionSpace = this->ghostMesh_[ghostMeshNo].get();
          std::array<int,2> ghostMeshCoordinates({neighbourCoordinatesLocal[0],0});
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
      else if (0 <= neighbourCoordinatesLocal[0] && neighbourCoordinatesLocal[0] <= this->meshPartition_->nElementsLocal(0)
        && 0 <= neighbourCoordinatesLocal[1] && neighbourCoordinatesLocal[1] <= this->meshPartition_->nElementsLocal(1))
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
    }  // x
  }   // y
  return false;
}

// 3D
template<typename MeshType, typename BasisFunctionType>
bool FunctionSpaceStructuredCheckNeighbouringElements<MeshType,BasisFunctionType,Mesh::isDim<3,MeshType>>::
checkNeighbouringElements(const Vec3 &point, element_no_t &elementNo, int &ghostMeshNo, std::array<double,MeshType::dim()> &xi)
{
  VLOG(1) << "checkNeighbouringElements(elementNo = " << elementNo << ", ghostMeshNo = " << ghostMeshNo << ", initial xi = " << xi;

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

  FunctionSpaceStructuredFindPositionBase<MeshType,BasisFunctionType> *functionSpace = this;

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
          if (this->ghostMesh_[(int)Mesh::face_t::face2Minus] == nullptr && this->ghostMesh_[(int)Mesh::face_t::face2Plus] == nullptr)
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
          if (this->ghostMesh_[(int)Mesh::face_t::face0Minus] != nullptr)
          {
            ghostMeshNo = (int)Mesh::face_t::face0Minus;
            functionSpace = this->ghostMesh_[ghostMeshNo].get();
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
          if (this->ghostMesh_[(int)Mesh::face_t::face0Plus] != nullptr)
          {
            ghostMeshNo = (int)Mesh::face_t::face0Plus;
            functionSpace = this->ghostMesh_[ghostMeshNo].get();
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
          if (this->ghostMesh_[(int)Mesh::face_t::face1Minus] != nullptr)
          {
            ghostMeshNo = (int)Mesh::face_t::face1Minus;
            functionSpace = this->ghostMesh_[ghostMeshNo].get();
            std::array<int,3> ghostMeshCoordinates({neighbourCoordinatesLocal[0],0,neighbourCoordinatesLocal[2]});
            neighbourElementNo = this->ghostMesh_[ghostMeshNo]->meshPartition()->getElementNoLocal(ghostMeshCoordinates);
            VLOG(1) << "1- neighbourElementNo: " << neighbourElementNo << ", ghostMeshCoordinates: " << ghostMeshCoordinates
              << " / (" << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(0) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(1) << "," << this->ghostMesh_[ghostMeshNo]->meshPartition()->nElementsLocal(2) << ")";
          }
        }
        else if (neighbourCoordinatesLocal[1] == this->meshPartition_->nElementsLocal(1))  // if at back boundary
        {
          if (this->ghostMesh_[(int)Mesh::face_t::face1Plus] != nullptr)
          {
            ghostMeshNo = (int)Mesh::face_t::face1Plus;
            functionSpace = this->ghostMesh_[ghostMeshNo].get();
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
          if (this->ghostMesh_[(int)Mesh::face_t::face2Minus] != nullptr)
          {
            ghostMeshNo = (int)Mesh::face_t::face2Minus;
            functionSpace = this->ghostMesh_[ghostMeshNo].get();
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
          if (this->ghostMesh_[(int)Mesh::face_t::face2Plus] != nullptr)
          {
            ghostMeshNo = (int)Mesh::face_t::face2Plus;
            functionSpace = this->ghostMesh_[ghostMeshNo].get();
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
      }  // x
    }  // y
  }   // z
  return false;
}

};  // namespace
