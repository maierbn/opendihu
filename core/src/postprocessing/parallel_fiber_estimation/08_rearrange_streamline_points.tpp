#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
rearrangeStreamlinePoints(std::vector<std::vector<Vec3>> &streamlineZPoints, std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &boundaryPointsSubdomain,
                          std::array<std::vector<Vec3>,4> &cornerStreamlines,
                          std::array<std::array<std::vector<bool>,4>,8> &boundaryPointsSubdomainAreValid,
                          std::array<bool,4> &subdomainIsAtBoundary)
{
  // boundaryPointsSubdomain[subdomain index][face_t][z-level][point index]

  // the numbering of the subdomains from 0-7 is as expected (morton numbering)

  // allocate space for boundaryPointsSubdomain
  for (int subdomainIndex = 0; subdomainIndex < 8; subdomainIndex++)
  {
    for (int faceNo = (int)Mesh::face_t::face0Minus; faceNo <= (int)Mesh::face_t::face1Plus; faceNo++)
    {
      boundaryPointsSubdomainAreValid[subdomainIndex][faceNo].resize(nBoundaryPointsX_,true);   // resize to number of points with same z level per face of subdomain
      boundaryPointsSubdomain[subdomainIndex][faceNo].resize(nBoundaryPointsZ_);   // resize to number of z levels
      for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZ_; zLevelIndex++)
      {
        boundaryPointsSubdomain[subdomainIndex][faceNo][zLevelIndex].resize(nBoundaryPointsX_,Vec3({0.0,0.0,0.0}));   // resize to number of points with same z level per face of subdomain
      }
    }
  }

  // assign sampled points to the data structure boundaryPointsSubdomain, which contains the points for each subdomain and face, as list of points for each z level

  int nBoundaryPointsXNew_ = nBoundaryPointsX_*2 - 1;
  int nBoundaryPointsZNew_ = nBoundaryPointsZ_*2 - 1;  // = (nBoundaryPointsZ_ - 1)*2 + 1
  int streamlineIndex = 0;

  // boundary indices for face0Minus and face0Plus (vertical)
  int iBeginVertical = 0;
  int iEndVertical = nBoundaryPointsXNew_;

  if (subdomainIsAtBoundary[(int)Mesh::face_t::face1Minus])
    iBeginVertical += 1;

  if (subdomainIsAtBoundary[(int)Mesh::face_t::face1Plus])
    iEndVertical -= 1;

  // boundary indices for face1Minus and face1Plus (horizontal)
  int iBeginHorizontal = 0;
  int iEndHorizontal = nBoundaryPointsXNew_;

  if (subdomainIsAtBoundary[(int)Mesh::face_t::face0Minus])
    iBeginHorizontal += 1;

  if (subdomainIsAtBoundary[(int)Mesh::face_t::face0Plus])
    iEndHorizontal -= 1;

  // face0Minus
  if (!subdomainIsAtBoundary[(int)Mesh::face_t::face0Minus])
  {
    // subdomains 0,4
    int pointIndex = iBeginVertical;
    for (int i = iBeginVertical; i < nBoundaryPointsX_; i++, streamlineIndex++, pointIndex++)
    {
      // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
      if (streamlineZPoints[streamlineIndex].size() == nBoundaryPointsZNew_)
      {
        // loop over bottom half of the streamline points
        for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZ_; zLevelIndex++)
        {
          const int subdomainIndex = 0;
          boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Minus][pointIndex] = true;
        }

        // loop over top half of the streamline points
        for (int zLevelIndex = nBoundaryPointsZ_-1; zLevelIndex < nBoundaryPointsZNew_; zLevelIndex++)
        {
          const int subdomainIndex = 4;
          const int zLevelIndexNew = zLevelIndex - (nBoundaryPointsZ_-1);
          boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Minus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Minus][pointIndex] = true;
        }
      }
      else
      {
        boundaryPointsSubdomainAreValid[0][(int)Mesh::face_t::face0Minus][pointIndex] = false;
        boundaryPointsSubdomainAreValid[4][(int)Mesh::face_t::face0Minus][pointIndex] = false;
        // set seed point at position 0
        boundaryPointsSubdomain[0][(int)Mesh::face_t::face0Minus][0][pointIndex] = streamlineZPoints[streamlineIndex][0];
      }
    }

    // subdomains 2,6
    pointIndex = 0;
    streamlineIndex--;  // consider center point twice
    for (int i = nBoundaryPointsX_-1; i < iEndVertical; i++, streamlineIndex++, pointIndex++)
    {
      // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
      if (streamlineZPoints[streamlineIndex].size() == nBoundaryPointsZNew_)
      {
        // loop over bottom half of the streamline points
        for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZ_; zLevelIndex++)
        {
          const int subdomainIndex = 2;
          boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Minus][pointIndex] = true;
        }

        // loop over top half of the streamline points
        for (int zLevelIndex = nBoundaryPointsZ_-1; zLevelIndex < nBoundaryPointsZNew_; zLevelIndex++)
        {
          const int subdomainIndex = 6;
          const int zLevelIndexNew = zLevelIndex - (nBoundaryPointsZ_-1);
          boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Minus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Minus][pointIndex] = true;
        }
      }
      else
      {
        boundaryPointsSubdomainAreValid[2][(int)Mesh::face_t::face0Minus][pointIndex] = false;
        boundaryPointsSubdomainAreValid[6][(int)Mesh::face_t::face0Minus][pointIndex] = false;
        // set seed point at position 0
        boundaryPointsSubdomain[2][(int)Mesh::face_t::face0Minus][0][pointIndex] = streamlineZPoints[streamlineIndex][0];
      }
    }
  }

  // face0Plus
  if (!subdomainIsAtBoundary[(int)Mesh::face_t::face0Plus])
  {
    // subdomains 1,5
    int pointIndex = iBeginVertical;
    for (int i = iBeginVertical; i < nBoundaryPointsX_; i++, streamlineIndex++, pointIndex++)
    {
      // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
      if (streamlineZPoints[streamlineIndex].size() == nBoundaryPointsZNew_)
      {
        // loop over bottom half of the streamline points
        for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZ_; zLevelIndex++)
        {
          const int subdomainIndex = 1;
          boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Plus][pointIndex] = true;
        }

        // loop over top half of the streamline points
        for (int zLevelIndex = nBoundaryPointsZ_-1; zLevelIndex < nBoundaryPointsZNew_; zLevelIndex++)
        {
          const int subdomainIndex = 5;
          const int zLevelIndexNew = zLevelIndex - (nBoundaryPointsZ_-1);
          boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Plus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Plus][pointIndex] = true;
        }
      }
      else
      {
        boundaryPointsSubdomainAreValid[1][(int)Mesh::face_t::face0Plus][pointIndex] = false;
        boundaryPointsSubdomainAreValid[5][(int)Mesh::face_t::face0Plus][pointIndex] = false;
        // set seed point at position 0
        boundaryPointsSubdomain[1][(int)Mesh::face_t::face0Plus][0][pointIndex] = streamlineZPoints[streamlineIndex][0];
      }
    }

    // subdomains 3,7
    pointIndex = 0;
    streamlineIndex--;  // consider center point twice
    for (int i = nBoundaryPointsX_-1; i < iEndVertical; i++, streamlineIndex++, pointIndex++)
    {
      // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
      if (streamlineZPoints[streamlineIndex].size() == nBoundaryPointsZNew_)
      {
        // loop over bottom half of the streamline points
        for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZ_; zLevelIndex++)
        {
          const int subdomainIndex = 3;
          boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Plus][pointIndex] = true;
        }

        // loop over top half of the streamline points
        for (int zLevelIndex = nBoundaryPointsZ_-1; zLevelIndex < nBoundaryPointsZNew_; zLevelIndex++)
        {
          const int subdomainIndex = 7;
          const int zLevelIndexNew = zLevelIndex - (nBoundaryPointsZ_-1);
          boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Plus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Plus][pointIndex] = true;
        }
      }
      else
      {
        boundaryPointsSubdomainAreValid[3][(int)Mesh::face_t::face0Plus][pointIndex] = false;
        boundaryPointsSubdomainAreValid[7][(int)Mesh::face_t::face0Plus][pointIndex] = false;
        // set seed point at position 0
        boundaryPointsSubdomain[3][(int)Mesh::face_t::face0Plus][0][pointIndex] = streamlineZPoints[streamlineIndex][0];
      }
    }
  }

  // face1Minus (without corner points)
  if (!subdomainIsAtBoundary[(int)Mesh::face_t::face1Minus])
  {
    // subdomains 0,4
    int pointIndex = iBeginHorizontal;
    for (int i = iBeginHorizontal; i < nBoundaryPointsX_; i++, streamlineIndex++, pointIndex++)
    {
      // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
      if (streamlineZPoints[streamlineIndex].size() == nBoundaryPointsZNew_)
      {
        // loop over bottom half of the streamline points
        for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZ_; zLevelIndex++)
        {
          const int subdomainIndex = 0;
          boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Minus][pointIndex] = true;
        }

        // loop over top half of the streamline points
        for (int zLevelIndex = nBoundaryPointsZ_-1; zLevelIndex < nBoundaryPointsZNew_; zLevelIndex++)
        {
          const int subdomainIndex = 4;
          const int zLevelIndexNew = zLevelIndex - (nBoundaryPointsZ_-1);
          boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Minus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Minus][pointIndex] = true;
        }
      }
      else
      {
        boundaryPointsSubdomainAreValid[0][(int)Mesh::face_t::face1Minus][pointIndex] = false;
        boundaryPointsSubdomainAreValid[4][(int)Mesh::face_t::face1Minus][pointIndex] = false;
        // set seed point at position 0
        boundaryPointsSubdomain[0][(int)Mesh::face_t::face1Minus][0][pointIndex] = streamlineZPoints[streamlineIndex][0];
      }
    }

    // subdomains 1,5
    pointIndex = 0;
    streamlineIndex--;  // consider center point twice
    for (int i = nBoundaryPointsX_-1; i < iEndHorizontal; i++, streamlineIndex++, pointIndex++)
    {
      // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
      if (streamlineZPoints[streamlineIndex].size() == nBoundaryPointsZNew_)
      {
        // loop over bottom half of the streamline points
        for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZ_; zLevelIndex++)
        {
          const int subdomainIndex = 1;
          boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Minus][pointIndex] = true;
        }

        // loop over top half of the streamline points
        for (int zLevelIndex = nBoundaryPointsZ_-1; zLevelIndex < nBoundaryPointsZNew_; zLevelIndex++)
        {
          const int subdomainIndex = 5;
          const int zLevelIndexNew = zLevelIndex - (nBoundaryPointsZ_-1);
          boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Minus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Minus][pointIndex] = true;
        }
      }
      else
      {
        boundaryPointsSubdomainAreValid[1][(int)Mesh::face_t::face1Minus][pointIndex] = false;
        boundaryPointsSubdomainAreValid[5][(int)Mesh::face_t::face1Minus][pointIndex] = false;
        // set seed point at position 0
        boundaryPointsSubdomain[1][(int)Mesh::face_t::face1Minus][0][pointIndex] = streamlineZPoints[streamlineIndex][0];
      }
    }
  }

  // face1Plus (with corner points)
  if (!subdomainIsAtBoundary[(int)Mesh::face_t::face1Plus])
  {
    // subdomains 2,6
    int pointIndex = iBeginHorizontal;
    for (int i = iBeginHorizontal; i < nBoundaryPointsX_; i++, streamlineIndex++, pointIndex++)
    {
      // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
      if (streamlineZPoints[streamlineIndex].size() == nBoundaryPointsZNew_)
      {
        // loop over bottom half of the streamline points
        for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZ_; zLevelIndex++)
        {
          const int subdomainIndex = 2;
          boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Plus][pointIndex] = true;
        }

        // loop over top half of the streamline points
        for (int zLevelIndex = nBoundaryPointsZ_-1; zLevelIndex < nBoundaryPointsZNew_; zLevelIndex++)
        {
          const int subdomainIndex = 6;
          const int zLevelIndexNew = zLevelIndex - (nBoundaryPointsZ_-1);
          boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Plus][pointIndex] = true;
        }
      }
      else
      {
        boundaryPointsSubdomainAreValid[2][(int)Mesh::face_t::face1Plus][pointIndex] = false;
        boundaryPointsSubdomainAreValid[6][(int)Mesh::face_t::face1Plus][pointIndex] = false;
        // set seed point at position 0
        boundaryPointsSubdomain[2][(int)Mesh::face_t::face1Plus][0][pointIndex] = streamlineZPoints[streamlineIndex][0];
      }
    }

    // subdomains 3,7
    pointIndex = 0;
    streamlineIndex--;  // consider center point twice
    for (int i = nBoundaryPointsX_-1; i < iEndHorizontal; i++, streamlineIndex++, pointIndex++)
    {
      // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
      if (streamlineZPoints[streamlineIndex].size() == nBoundaryPointsZNew_)
      {
        // loop over bottom half of the streamline points
        for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZ_; zLevelIndex++)
        {
          const int subdomainIndex = 3;
          boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Plus][pointIndex] = true;
        }

        // loop over top half of the streamline points
        for (int zLevelIndex = nBoundaryPointsZ_-1; zLevelIndex < nBoundaryPointsZNew_; zLevelIndex++)
        {
          const int subdomainIndex = 7;
          const int zLevelIndexNew = zLevelIndex - (nBoundaryPointsZ_-1);
          boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Plus][pointIndex] = true;
        }
      }
      else
      {
        boundaryPointsSubdomainAreValid[3][(int)Mesh::face_t::face1Plus][pointIndex] = false;
        boundaryPointsSubdomainAreValid[7][(int)Mesh::face_t::face1Plus][pointIndex] = false;
        // set seed point at position 0
        boundaryPointsSubdomain[3][(int)Mesh::face_t::face1Plus][0][pointIndex] = streamlineZPoints[streamlineIndex][0];
      }
    }
  }

  LOG(DEBUG) << "starting with horizontal center line, streamlineIndex = " << streamlineIndex;
  // horizontal center line (with corner points)
  // subdomains 0,4
  int pointIndex = iBeginHorizontal;
  int streamlineIndexHorizontalStart = streamlineIndex;
  for (int i = iBeginHorizontal; i < nBoundaryPointsX_; i++, streamlineIndex++, pointIndex++)
  {
    // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
    if (streamlineZPoints[streamlineIndex].size() == nBoundaryPointsZNew_)
    {
      // loop over bottom half of the streamline points
      for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZ_; zLevelIndex++)
      {
        const int subdomainIndex = 0;
        VLOG(1) << boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus].size() << " levels ";
        VLOG(1) << boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndex].size() << " points";
        VLOG(1) << "z: " << zLevelIndex << ", i=" << i;
        boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Plus][pointIndex] = true;
      }

      // loop over top half of the streamline points
      for (int zLevelIndex = nBoundaryPointsZ_-1; zLevelIndex < nBoundaryPointsZNew_; zLevelIndex++)
      {
        const int subdomainIndex = 4;
        const int zLevelIndexNew = zLevelIndex - (nBoundaryPointsZ_-1);
        boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Plus][pointIndex] = true;
      }
    }
    else
    {
      boundaryPointsSubdomainAreValid[0][(int)Mesh::face_t::face1Plus][pointIndex] = false;
      boundaryPointsSubdomainAreValid[4][(int)Mesh::face_t::face1Plus][pointIndex] = false;
      // set seed point at position 0
      boundaryPointsSubdomain[0][(int)Mesh::face_t::face1Plus][0][pointIndex] = streamlineZPoints[streamlineIndex][0];
    }
  }

  // subdomains 1,5
  pointIndex = 0;
  streamlineIndex--;  // consider center point twice
  for (int i = nBoundaryPointsX_-1; i < iEndHorizontal; i++, streamlineIndex++, pointIndex++)
  {
    // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
    if (streamlineZPoints[streamlineIndex].size() == nBoundaryPointsZNew_)
    {
      // loop over bottom half of the streamline points
      for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZ_; zLevelIndex++)
      {
        const int subdomainIndex = 1;
        boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Plus][pointIndex] = true;
      }

      // loop over top half of the streamline points
      for (int zLevelIndex = nBoundaryPointsZ_-1; zLevelIndex < nBoundaryPointsZNew_; zLevelIndex++)
      {
        const int subdomainIndex = 5;
        const int zLevelIndexNew = zLevelIndex - (nBoundaryPointsZ_-1);
        boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Plus][pointIndex] = true;
      }
    }
    else
    {
      boundaryPointsSubdomainAreValid[1][(int)Mesh::face_t::face1Plus][pointIndex] = false;
      boundaryPointsSubdomainAreValid[5][(int)Mesh::face_t::face1Plus][pointIndex] = false;
      // set seed point at position 0
      boundaryPointsSubdomain[1][(int)Mesh::face_t::face1Plus][0][pointIndex] = streamlineZPoints[streamlineIndex][0];
    }
  }

  // subdomains 2,6
  pointIndex = iBeginHorizontal;
  streamlineIndex = streamlineIndexHorizontalStart;
  for (int i = iBeginHorizontal; i < nBoundaryPointsX_; i++, streamlineIndex++, pointIndex++)
  {
    // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
    if (streamlineZPoints[streamlineIndex].size() == nBoundaryPointsZNew_)
    {
      // loop over bottom half of the streamline points
      for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZ_; zLevelIndex++)
      {
        const int subdomainIndex = 2;
        boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Minus][pointIndex] = true;
      }

      // loop over top half of the streamline points
      for (int zLevelIndex = nBoundaryPointsZ_-1; zLevelIndex < nBoundaryPointsZNew_; zLevelIndex++)
      {
        const int subdomainIndex = 6;
        const int zLevelIndexNew = zLevelIndex - (nBoundaryPointsZ_-1);
        boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Minus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Minus][pointIndex] = true;
      }
    }
    else
    {
      boundaryPointsSubdomainAreValid[2][(int)Mesh::face_t::face1Minus][pointIndex] = false;
      boundaryPointsSubdomainAreValid[6][(int)Mesh::face_t::face1Minus][pointIndex] = false;
      // set seed point at position 0
      boundaryPointsSubdomain[2][(int)Mesh::face_t::face1Minus][0][pointIndex] = streamlineZPoints[streamlineIndex][0];
    }
  }

  // subdomains 3,7
  pointIndex = 0;
  streamlineIndex--;  // consider center point twice
  for (int i = nBoundaryPointsX_-1; i < iEndHorizontal; i++, streamlineIndex++, pointIndex++)
  {
    // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
    if (streamlineZPoints[streamlineIndex].size() == nBoundaryPointsZNew_)
    {

      // loop over bottom half of the streamline points
      for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZ_; zLevelIndex++)
      {
        const int subdomainIndex = 3;
        boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Minus][pointIndex] = true;
      }

      // loop over top half of the streamline points
      for (int zLevelIndex = nBoundaryPointsZ_-1; zLevelIndex < nBoundaryPointsZNew_; zLevelIndex++)
      {
        const int subdomainIndex = 7;
        const int zLevelIndexNew = zLevelIndex - (nBoundaryPointsZ_-1);
        boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Minus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Minus][pointIndex] = true;
      }
    }
    else
    {
      boundaryPointsSubdomainAreValid[3][(int)Mesh::face_t::face1Minus][pointIndex] = false;
      boundaryPointsSubdomainAreValid[7][(int)Mesh::face_t::face1Minus][pointIndex] = false;
      // set seed point at position 0
      boundaryPointsSubdomain[3][(int)Mesh::face_t::face1Minus][0][pointIndex] = streamlineZPoints[streamlineIndex][0];
    }
  }

  LOG(DEBUG) << "starting with vertical center line, streamlineIndex = " << streamlineIndex;

  // vertical center line (with corner points and center point)
  // subdomains 0,4, 1,5
  pointIndex = iBeginVertical;
  for (int i = iBeginVertical; i < nBoundaryPointsX_; i++, streamlineIndex++, pointIndex++)
  {
    // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
    if (streamlineZPoints[streamlineIndex].size() == nBoundaryPointsZNew_)
    {
      // subdomains 0,4,
      // loop over bottom half of the streamline points
      for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZ_; zLevelIndex++)
      {
        const int subdomainIndex = 0;
        boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Plus][pointIndex] = true;
      }

      // loop over top half of the streamline points
      for (int zLevelIndex = nBoundaryPointsZ_-1; zLevelIndex < nBoundaryPointsZNew_; zLevelIndex++)
      {
        const int subdomainIndex = 4;
        const int zLevelIndexNew = zLevelIndex - (nBoundaryPointsZ_-1);
        boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Plus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Plus][pointIndex] = true;
      }
    }
    else
    {
      boundaryPointsSubdomainAreValid[0][(int)Mesh::face_t::face0Plus][pointIndex] = false;
      boundaryPointsSubdomainAreValid[4][(int)Mesh::face_t::face0Plus][pointIndex] = false;
      // set seed point at position 0
      boundaryPointsSubdomain[0][(int)Mesh::face_t::face0Plus][0][pointIndex] = streamlineZPoints[streamlineIndex][0];
    }

    // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
    if (streamlineZPoints[streamlineIndex].size() == nBoundaryPointsZNew_)
    {
      // subdomains 1,5
      // loop over bottom half of the streamline points
      for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZ_; zLevelIndex++)
      {
        const int subdomainIndex = 1;
        boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Minus][pointIndex] = true;
      }

      // loop over top half of the streamline points
      for (int zLevelIndex = nBoundaryPointsZ_-1; zLevelIndex < nBoundaryPointsZNew_; zLevelIndex++)
      {
        const int subdomainIndex = 5;
        const int zLevelIndexNew = zLevelIndex - (nBoundaryPointsZ_-1);
        boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Minus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Minus][pointIndex] = true;
      }
    }
    else
    {
      boundaryPointsSubdomainAreValid[1][(int)Mesh::face_t::face0Minus][pointIndex] = false;
      boundaryPointsSubdomainAreValid[5][(int)Mesh::face_t::face0Minus][pointIndex] = false;
      // set seed point at position 0
      boundaryPointsSubdomain[1][(int)Mesh::face_t::face0Minus][0][pointIndex] = streamlineZPoints[streamlineIndex][0];
    }
  }

  LOG(DEBUG) << "nBoundaryPointsZ_: " << nBoundaryPointsZ_ << ", nBoundaryPointsX_: " << nBoundaryPointsX_ << ", iEndVertical: " << iEndVertical;

  // subdomains 2,6, 3,7
  pointIndex = 0;
  streamlineIndex--;  // consider center point twice
  for (int i = nBoundaryPointsX_-1; i < iEndVertical; i++, streamlineIndex++, pointIndex++)
  {
    // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
    if (streamlineZPoints[streamlineIndex].size() == nBoundaryPointsZNew_)
    {
      VLOG(1) << "";
      VLOG(1) << "i=" << i << " (iEndVertical=" << iEndVertical << ")";
      VLOG(1) << "streamlineIndex: " << streamlineIndex;
      VLOG(1) << "pointIndex: " << pointIndex << ", nBoundaryPointsZ_: " << nBoundaryPointsZ_;

      // subdomains 2,6
      // loop over bottom half of the streamline points
      for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZ_; zLevelIndex++)
      {
        const int subdomainIndex = 2;
        VLOG(1) << "i=" << i << " (iEndVertical=" << iEndVertical << "), set boundaryPointsSubdomain[" << subdomainIndex << "]["
          << (int)Mesh::face_t::face0Plus << "][" << zLevelIndex << "][" <<pointIndex << "] = streamlineZPoints[" << streamlineIndex
          << "][" << zLevelIndex << "] = " << streamlineZPoints[streamlineIndex][zLevelIndex];
        boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Plus][pointIndex] = true;
      }

      // loop over top half of the streamline points
      for (int zLevelIndex = nBoundaryPointsZ_-1; zLevelIndex < nBoundaryPointsZNew_; zLevelIndex++)
      {
        const int subdomainIndex = 6;
        const int zLevelIndexNew = zLevelIndex - (nBoundaryPointsZ_-1);
        boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Plus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Plus][pointIndex] = true;
      }
    }
    else
    {
      boundaryPointsSubdomainAreValid[2][(int)Mesh::face_t::face0Plus][pointIndex] = false;
      boundaryPointsSubdomainAreValid[6][(int)Mesh::face_t::face0Plus][pointIndex] = false;
      // set seed point at position 0
      boundaryPointsSubdomain[2][(int)Mesh::face_t::face0Plus][0][pointIndex] = streamlineZPoints[streamlineIndex][0];
    }

    // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
    if (streamlineZPoints[streamlineIndex].size() == nBoundaryPointsZNew_)
    {
      // subdomains 3,7
      // loop over bottom half of the streamline points
      for (int zLevelIndex = 0; zLevelIndex < nBoundaryPointsZ_; zLevelIndex++)
      {
        const int subdomainIndex = 3;
        boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Minus][pointIndex] = true;
      }

      // loop over top half of the streamline points
      for (int zLevelIndex = nBoundaryPointsZ_-1; zLevelIndex < nBoundaryPointsZNew_; zLevelIndex++)
      {
        const int subdomainIndex = 7;
        const int zLevelIndexNew = zLevelIndex - (nBoundaryPointsZ_-1);
        boundaryPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Minus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        boundaryPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Minus][pointIndex] = true;
      }
    }
    else
    {
      boundaryPointsSubdomainAreValid[3][(int)Mesh::face_t::face0Minus][pointIndex] = false;
      boundaryPointsSubdomainAreValid[7][(int)Mesh::face_t::face0Minus][pointIndex] = false;
      // set seed point at position 0
      boundaryPointsSubdomain[3][(int)Mesh::face_t::face0Minus][0][pointIndex] = streamlineZPoints[streamlineIndex][0];
    }
  }

  for (int cornerStreamlineIndex = 0; cornerStreamlineIndex != 4; cornerStreamlineIndex++, streamlineIndex++)
  {
    cornerStreamlines[cornerStreamlineIndex].resize(streamlineZPoints[streamlineIndex].size());
    for (int zLevelIndex = 0; zLevelIndex < streamlineZPoints[streamlineIndex].size(); zLevelIndex++)
    {
      cornerStreamlines[cornerStreamlineIndex][zLevelIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
    }
  }

#ifndef NDEBUG
  std::stringstream s;
  for (int subdomainIndex = 0; subdomainIndex < 8; subdomainIndex++)
  {
    s.str("");
    s << "06_boundary_points_subdomain_" << subdomainIndex;
    PyObject_CallFunction(functionOutputBoundaryPoints_, "s i i O f", s.str().c_str(), currentRankSubset_->ownRankNo(), level_,
                          PythonUtility::convertToPython<std::array<std::vector<std::vector<Vec3>>,4>>::get(boundaryPointsSubdomain[subdomainIndex]), 0.1);
    PythonUtility::checkForError();
  }
  PyObject_CallFunction(functionOutputStreamlines_, "s i i O f", "06_corner_streamlines", currentRankSubset_->ownRankNo(), level_,
                        PythonUtility::convertToPython<std::array<std::vector<Vec3>,4>>::get(cornerStreamlines), 0.1);
  PythonUtility::checkForError();
#endif
}

} // namespace
