#include "function_space/09_function_space_find_position_base.h"

#include <cmath>
#include <array>
#include <sstream>

#include "easylogging++.h"
#include "mesh/face_t.h"

namespace FunctionSpace
{

// composite mesh
template<int D,typename BasisFunctionType>
bool FunctionSpaceStructuredFindPositionBase<Mesh::CompositeOfDimension<D>,BasisFunctionType>::
findPosition(Vec3 point, element_no_t &elementNoLocal, int &ghostMeshNo, std::array<double,MeshType::dim()> &xi, bool startSearchInCurrentElement, double xiTolerance)
{
  const element_no_t nElements = this->nElementsLocal();
  VLOG(2) << "findPosition, elementNoLocal: " << elementNoLocal << ", ghostMeshNo: " << ghostMeshNo
    << ", startSearchInCurrentElement: " << startSearchInCurrentElement << ", xiTolerance: " << xiTolerance << ", xi: " << xi;

  // set starting no to 0 if it was not given and is thus arbitrarily initialized
  if (elementNoLocal < 0 || elementNoLocal >= nElements)
    elementNoLocal = 0;

  int subMeshOfElementSubmeshNo = 0;

  // if the search should start in the current element given by elementNoLocal
  if (startSearchInCurrentElement)
  {
    FunctionSpaceStructuredFindPositionBase<MeshType,BasisFunctionType> *functionSpace = this;

    // check if point is already in current element
    if (functionSpace->pointIsInElement(point, elementNoLocal, xi, xiTolerance))
    {

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

      return true;
    }

    // point is not in current element, consider the neighbouring elements and ghost meshes

    VLOG(2) << "point is not in current element, now check neighbouring elements";

    // search in the current submesh
    element_no_t elementOnMeshNoLocal = 0;
    this->meshPartition_->getSubMeshNoAndElementNoLocal(elementNoLocal, subMeshOfElementSubmeshNo, elementOnMeshNoLocal);

    bool nodeFound = this->subFunctionSpaces_[subMeshOfElementSubmeshNo]->findPosition(point, elementOnMeshNoLocal, ghostMeshNo, xi,
                                                                                       startSearchInCurrentElement, xiTolerance);
    if (nodeFound)
    {
      return true;
    }
    else
    {
      VLOG(2) << "(returning from recursion), Point was not found in the entire submesh no " << subMeshOfElementSubmeshNo << ", xi: " << xi << ", now check all elements";
    }
  }

  // search among all elements in submeshes other than subMeshNo (because it was already searched there)
  for (int subMeshNo = 0; subMeshNo < this->subFunctionSpaces_.size(); subMeshNo++)
  {
    if (subMeshNo == subMeshOfElementSubmeshNo)
      continue;

    bool nodeFound = this->subFunctionSpaces_[subMeshNo]->findPosition(point, elementNoLocal, ghostMeshNo, xi,
                                                                       false, xiTolerance);
    if (nodeFound)
    {
      return true;
    }
  }

  VLOG(1) << "Could not find any containing element";
  return false;
}

} // namespace
