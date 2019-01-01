#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
createSeedPoints(const std::array<bool,4> &subdomainIsAtBorder, int seedPointsZIndex, const std::vector<Vec3> &nodePositions, std::vector<Vec3> &seedPoints)
{
  int subdomainNNodesX = nBorderPointsXNew_;
  int subdomainNNodesY = nBorderPointsXNew_;

  // boundary indices for face0Minus and face0Plus (vertical direction)
  int iBeginVertical = 0;
  int iEndVertical = nBorderPointsXNew_;

  if (subdomainIsAtBorder[(int)Mesh::face_t::face1Minus])
    iBeginVertical += 1;

  if (subdomainIsAtBorder[(int)Mesh::face_t::face1Plus])
    iEndVertical -= 1;

  // boundary indices for face1Minus and face1Plus (horizontal direction)
  int iBeginHorizontal = 0;
  int iEndHorizontal = nBorderPointsXNew_;

  if (subdomainIsAtBorder[(int)Mesh::face_t::face0Minus])
    iBeginHorizontal += 1;

  if (subdomainIsAtBorder[(int)Mesh::face_t::face0Plus])
    iEndHorizontal -= 1;

  // face0Minus
  if (!subdomainIsAtBorder[(int)Mesh::face_t::face0Minus])
  {
    for (int i = iBeginVertical; i < iEndVertical; i++)
    {
      seedPoints.push_back(nodePositions[seedPointsZIndex*subdomainNNodesX*subdomainNNodesY + i*subdomainNNodesX + 0]);
    }
  }

  // face0Plus
  if (!subdomainIsAtBorder[(int)Mesh::face_t::face0Plus])
  {
    for (int i = iBeginVertical; i < iEndVertical; i++)
    {
      seedPoints.push_back(nodePositions[seedPointsZIndex*subdomainNNodesX*subdomainNNodesY + i*subdomainNNodesX + (subdomainNNodesX-1)]);
    }
  }

  // face1Minus (with corner points)
  if (!subdomainIsAtBorder[(int)Mesh::face_t::face1Minus])
  {
    for (int i = iBeginHorizontal; i < iEndHorizontal; i++)
    {
      seedPoints.push_back(nodePositions[seedPointsZIndex*subdomainNNodesX*subdomainNNodesY + i]);
    }
  }

  // face1Plus (with corner points)
  if (!subdomainIsAtBorder[(int)Mesh::face_t::face1Plus])
  {
    for (int i = iBeginHorizontal; i < iEndHorizontal; i++)
    {
      seedPoints.push_back(nodePositions[seedPointsZIndex*subdomainNNodesX*subdomainNNodesY + (subdomainNNodesY-1)*subdomainNNodesX + i]);
    }
  }


  LOG(DEBUG) << "seedPoints: starting with horizontal center line, streamlineIndex = " << seedPoints.size();

  // horizontal center line (with corner points)
  for (int i = iBeginHorizontal; i < iEndHorizontal; i++)
  {
    seedPoints.push_back(nodePositions[seedPointsZIndex*subdomainNNodesX*subdomainNNodesY + int(subdomainNNodesY/2)*subdomainNNodesX + i]);
  }

  LOG(DEBUG) << "seedPoints: starting with vertical center line, streamlineIndex = " << seedPoints.size();

  // vertical center line (with corner points and center point)
  for (int i = iBeginVertical; i < iEndVertical; i++)
  {
    seedPoints.push_back(nodePositions[seedPointsZIndex*subdomainNNodesX*subdomainNNodesY + i*subdomainNNodesX + int(subdomainNNodesX/2)]);
  }

  LOG(DEBUG) << "seedPoints: end, streamlineIndex = " << seedPoints.size();

}

};  // namespace
