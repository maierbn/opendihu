#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
rearrangeStreamlinePoints(std::vector<std::vector<Vec3>> &streamlineZPoints, std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &borderPointsSubdomain,
                          std::array<std::array<std::vector<bool>,4>,8> &borderPointsSubdomainAreValid,
                          std::array<bool,4> &subdomainIsAtBorder)
{
  // borderPointsSubdomain[subdomain index][face_t][z-level][point index]

  // the numbering of the subdomains from 0-7 is as expected (morton numbering)

  // allocate space for borderPointsSubdomain
  for (int subdomainIndex = 0; subdomainIndex < 8; subdomainIndex++)
  {
    for (int faceNo = (int)Mesh::face_t::face0Minus; faceNo <= (int)Mesh::face_t::face1Plus; faceNo++)
    {
      borderPointsSubdomainAreValid[subdomainIndex][faceNo].resize(nBorderPointsX_,true);   // resize to number of points with same z level per face of subdomain
      borderPointsSubdomain[subdomainIndex][faceNo].resize(nBorderPointsZ_);   // resize to number of z levels
      for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
      {
        borderPointsSubdomain[subdomainIndex][faceNo][zLevelIndex].resize(nBorderPointsX_);   // resize to number of points with same z level per face of subdomain
      }
    }
  }

  // assign sampled points to the data structure borderPointsSubdomain, which contains the points for each subdomain and face, as list of points for each z level

  int nBorderPointsXNew_ = nBorderPointsX_*2 - 1;
  int nBorderPointsZNew_ = nBorderPointsZ_*2 - 1;  // = (nBorderPointsZ_ - 1)*2 + 1
  int streamlineIndex = 0;

  // boundary indices for face0Minus and face0Plus (vertical)
  int iBeginVertical = 0;
  int iEndVertical = nBorderPointsXNew_;

  if (subdomainIsAtBorder[(int)Mesh::face_t::face1Minus])
    iBeginVertical += 1;

  if (subdomainIsAtBorder[(int)Mesh::face_t::face1Plus])
    iEndVertical -= 1;

  // boundary indices for face1Minus and face1Plus (horizontal)
  int iBeginHorizontal = 0;
  int iEndHorizontal = nBorderPointsXNew_;

  if (subdomainIsAtBorder[(int)Mesh::face_t::face0Minus])
    iBeginHorizontal += 1;

  if (subdomainIsAtBorder[(int)Mesh::face_t::face0Plus])
    iEndHorizontal -= 1;

  // face0Minus
  if (!subdomainIsAtBorder[(int)Mesh::face_t::face0Minus])
  {
    // subdomains 0,4
    int pointIndex = iBeginVertical;
    for (int i = iBeginVertical; i < nBorderPointsX_; i++, streamlineIndex++, pointIndex++)
    {
      // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
      if (streamlineZPoints[streamlineIndex].size() == nBorderPointsZNew_)
      {
        // loop over bottom half of the streamline points
        for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
        {
          const int subdomainIndex = 0;
          borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Minus][pointIndex] = true;
        }

        // loop over top half of the streamline points
        for (int zLevelIndex = nBorderPointsZ_-1; zLevelIndex < nBorderPointsZNew_; zLevelIndex++)
        {
          const int subdomainIndex = 4;
          const int zLevelIndexNew = zLevelIndex - (nBorderPointsZ_-1);
          borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Minus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Minus][pointIndex] = true;
        }
      }
      else
      {
        borderPointsSubdomainAreValid[0][(int)Mesh::face_t::face0Minus][pointIndex] = false;
        borderPointsSubdomainAreValid[4][(int)Mesh::face_t::face0Minus][pointIndex] = false;
      }
    }

    // subdomains 2,6
    pointIndex = 0;
    streamlineIndex--;  // consider center point twice
    for (int i = nBorderPointsX_-1; i < iEndVertical; i++, streamlineIndex++, pointIndex++)
    {
      // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
      if (streamlineZPoints[streamlineIndex].size() == nBorderPointsZNew_)
      {
        // loop over bottom half of the streamline points
        for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
        {
          const int subdomainIndex = 2;
          borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Minus][pointIndex] = true;
        }

        // loop over top half of the streamline points
        for (int zLevelIndex = nBorderPointsZ_-1; zLevelIndex < nBorderPointsZNew_; zLevelIndex++)
        {
          const int subdomainIndex = 6;
          const int zLevelIndexNew = zLevelIndex - (nBorderPointsZ_-1);
          borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Minus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Minus][pointIndex] = true;
        }
      }
      else
      {
        borderPointsSubdomainAreValid[2][(int)Mesh::face_t::face0Minus][pointIndex] = false;
        borderPointsSubdomainAreValid[6][(int)Mesh::face_t::face0Minus][pointIndex] = false;
      }
    }
  }

  // face0Plus
  if (!subdomainIsAtBorder[(int)Mesh::face_t::face0Plus])
  {
    // subdomains 1,5
    int pointIndex = iBeginVertical;
    for (int i = iBeginVertical; i < nBorderPointsX_; i++, streamlineIndex++, pointIndex++)
    {
      // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
      if (streamlineZPoints[streamlineIndex].size() == nBorderPointsZNew_)
      {
        // loop over bottom half of the streamline points
        for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
        {
          const int subdomainIndex = 1;
          borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Plus][pointIndex] = true;
        }

        // loop over top half of the streamline points
        for (int zLevelIndex = nBorderPointsZ_-1; zLevelIndex < nBorderPointsZNew_; zLevelIndex++)
        {
          const int subdomainIndex = 5;
          const int zLevelIndexNew = zLevelIndex - (nBorderPointsZ_-1);
          borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Plus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Plus][pointIndex] = true;
        }
      }
      else
      {
        borderPointsSubdomainAreValid[1][(int)Mesh::face_t::face0Plus][pointIndex] = false;
        borderPointsSubdomainAreValid[5][(int)Mesh::face_t::face0Plus][pointIndex] = false;
      }
    }

    // subdomains 3,7
    pointIndex = 0;
    streamlineIndex--;  // consider center point twice
    for (int i = nBorderPointsX_-1; i < iEndVertical; i++, streamlineIndex++, pointIndex++)
    {
      // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
      if (streamlineZPoints[streamlineIndex].size() == nBorderPointsZNew_)
      {
        // loop over bottom half of the streamline points
        for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
        {
          const int subdomainIndex = 3;
          borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Plus][pointIndex] = true;
        }

        // loop over top half of the streamline points
        for (int zLevelIndex = nBorderPointsZ_-1; zLevelIndex < nBorderPointsZNew_; zLevelIndex++)
        {
          const int subdomainIndex = 7;
          const int zLevelIndexNew = zLevelIndex - (nBorderPointsZ_-1);
          borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Plus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Plus][pointIndex] = true;
        }
      }
      else
      {
        borderPointsSubdomainAreValid[3][(int)Mesh::face_t::face0Plus][pointIndex] = false;
        borderPointsSubdomainAreValid[7][(int)Mesh::face_t::face0Plus][pointIndex] = false;
      }
    }
  }

  // face1Minus (without corner points)
  if (!subdomainIsAtBorder[(int)Mesh::face_t::face1Minus])
  {
    // subdomains 0,4
    int pointIndex = iBeginHorizontal;
    for (int i = iBeginHorizontal; i < nBorderPointsX_; i++, streamlineIndex++, pointIndex++)
    {
      // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
      if (streamlineZPoints[streamlineIndex].size() == nBorderPointsZNew_)
      {
        // loop over bottom half of the streamline points
        for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
        {
          const int subdomainIndex = 0;
          borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Minus][pointIndex] = true;
        }

        // loop over top half of the streamline points
        for (int zLevelIndex = nBorderPointsZ_-1; zLevelIndex < nBorderPointsZNew_; zLevelIndex++)
        {
          const int subdomainIndex = 4;
          const int zLevelIndexNew = zLevelIndex - (nBorderPointsZ_-1);
          borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Minus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Minus][pointIndex] = true;
        }
      }
      else
      {
        borderPointsSubdomainAreValid[0][(int)Mesh::face_t::face1Minus][pointIndex] = false;
        borderPointsSubdomainAreValid[4][(int)Mesh::face_t::face1Minus][pointIndex] = false;
      }
    }

    // subdomains 1,5
    pointIndex = 0;
    streamlineIndex--;  // consider center point twice
    for (int i = nBorderPointsX_-1; i < iEndHorizontal; i++, streamlineIndex++, pointIndex++)
    {
      // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
      if (streamlineZPoints[streamlineIndex].size() == nBorderPointsZNew_)
      {
        // loop over bottom half of the streamline points
        for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
        {
          const int subdomainIndex = 1;
          borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Minus][pointIndex] = true;
        }

        // loop over top half of the streamline points
        for (int zLevelIndex = nBorderPointsZ_-1; zLevelIndex < nBorderPointsZNew_; zLevelIndex++)
        {
          const int subdomainIndex = 5;
          const int zLevelIndexNew = zLevelIndex - (nBorderPointsZ_-1);
          borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Minus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Minus][pointIndex] = true;
        }
      }
      else
      {
        borderPointsSubdomainAreValid[1][(int)Mesh::face_t::face1Minus][pointIndex] = false;
        borderPointsSubdomainAreValid[5][(int)Mesh::face_t::face1Minus][pointIndex] = false;
      }
    }
  }

  // face1Plus (with corner points)
  if (!subdomainIsAtBorder[(int)Mesh::face_t::face1Plus])
  {
    // subdomains 2,6
    int pointIndex = iBeginHorizontal;
    for (int i = iBeginHorizontal; i < nBorderPointsX_; i++, streamlineIndex++, pointIndex++)
    {
      // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
      if (streamlineZPoints[streamlineIndex].size() == nBorderPointsZNew_)
      {
        // loop over bottom half of the streamline points
        for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
        {
          const int subdomainIndex = 2;
          borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Plus][pointIndex] = true;
        }

        // loop over top half of the streamline points
        for (int zLevelIndex = nBorderPointsZ_-1; zLevelIndex < nBorderPointsZNew_; zLevelIndex++)
        {
          const int subdomainIndex = 6;
          const int zLevelIndexNew = zLevelIndex - (nBorderPointsZ_-1);
          borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Plus][pointIndex] = true;
        }
      }
      else
      {
        borderPointsSubdomainAreValid[2][(int)Mesh::face_t::face1Plus][pointIndex] = false;
        borderPointsSubdomainAreValid[6][(int)Mesh::face_t::face1Plus][pointIndex] = false;
      }
    }

    // subdomains 3,7
    pointIndex = 0;
    streamlineIndex--;  // consider center point twice
    for (int i = nBorderPointsX_-1; i < iEndHorizontal; i++, streamlineIndex++, pointIndex++)
    {
      // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
      if (streamlineZPoints[streamlineIndex].size() == nBorderPointsZNew_)
      {
        // loop over bottom half of the streamline points
        for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
        {
          const int subdomainIndex = 3;
          borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Plus][pointIndex] = true;
        }

        // loop over top half of the streamline points
        for (int zLevelIndex = nBorderPointsZ_-1; zLevelIndex < nBorderPointsZNew_; zLevelIndex++)
        {
          const int subdomainIndex = 7;
          const int zLevelIndexNew = zLevelIndex - (nBorderPointsZ_-1);
          borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
          borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Plus][pointIndex] = true;
        }
      }
      else
      {
        borderPointsSubdomainAreValid[3][(int)Mesh::face_t::face1Plus][pointIndex] = false;
        borderPointsSubdomainAreValid[7][(int)Mesh::face_t::face1Plus][pointIndex] = false;
      }
    }
  }

  LOG(DEBUG) << "starting with horizontal center line, streamlineIndex = " << streamlineIndex;
  // horizontal center line (with corner points)
  // subdomains 0,4
  int pointIndex = iBeginHorizontal;
  int streamlineIndexHorizontalStart = streamlineIndex;
  for (int i = iBeginHorizontal; i < nBorderPointsX_; i++, streamlineIndex++, pointIndex++)
  {
    // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
    if (streamlineZPoints[streamlineIndex].size() == nBorderPointsZNew_)
    {
      // loop over bottom half of the streamline points
      for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
      {
        const int subdomainIndex = 0;
        VLOG(1) << borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus].size() << " levels ";
        VLOG(1) << borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndex].size() << " points";
        VLOG(1) << "z: " << zLevelIndex << ", i=" << i;
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Plus][pointIndex] = true;
      }

      // loop over top half of the streamline points
      for (int zLevelIndex = nBorderPointsZ_-1; zLevelIndex < nBorderPointsZNew_; zLevelIndex++)
      {
        const int subdomainIndex = 4;
        const int zLevelIndexNew = zLevelIndex - (nBorderPointsZ_-1);
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Plus][pointIndex] = true;
      }
    }
    else
    {
      borderPointsSubdomainAreValid[0][(int)Mesh::face_t::face1Plus][pointIndex] = false;
      borderPointsSubdomainAreValid[4][(int)Mesh::face_t::face1Plus][pointIndex] = false;
    }
  }

  // subdomains 1,5
  pointIndex = 0;
  streamlineIndex--;  // consider center point twice
  for (int i = nBorderPointsX_-1; i < iEndHorizontal; i++, streamlineIndex++, pointIndex++)
  {
    // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
    if (streamlineZPoints[streamlineIndex].size() == nBorderPointsZNew_)
    {
      // loop over bottom half of the streamline points
      for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
      {
        const int subdomainIndex = 1;
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Plus][pointIndex] = true;
      }

      // loop over top half of the streamline points
      for (int zLevelIndex = nBorderPointsZ_-1; zLevelIndex < nBorderPointsZNew_; zLevelIndex++)
      {
        const int subdomainIndex = 5;
        const int zLevelIndexNew = zLevelIndex - (nBorderPointsZ_-1);
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Plus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Plus][pointIndex] = true;
      }
    }
    else
    {
      borderPointsSubdomainAreValid[1][(int)Mesh::face_t::face1Plus][pointIndex] = false;
      borderPointsSubdomainAreValid[5][(int)Mesh::face_t::face1Plus][pointIndex] = false;
    }
  }

  // subdomains 2,6
  pointIndex = iBeginHorizontal;
  streamlineIndex = streamlineIndexHorizontalStart;
  for (int i = iBeginHorizontal; i < nBorderPointsX_; i++, streamlineIndex++, pointIndex++)
  {
    // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
    if (streamlineZPoints[streamlineIndex].size() == nBorderPointsZNew_)
    {
      // loop over bottom half of the streamline points
      for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
      {
        const int subdomainIndex = 2;
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Minus][pointIndex] = true;
      }

      // loop over top half of the streamline points
      for (int zLevelIndex = nBorderPointsZ_-1; zLevelIndex < nBorderPointsZNew_; zLevelIndex++)
      {
        const int subdomainIndex = 6;
        const int zLevelIndexNew = zLevelIndex - (nBorderPointsZ_-1);
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Minus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Minus][pointIndex] = true;
      }
    }
    else
    {
      borderPointsSubdomainAreValid[2][(int)Mesh::face_t::face1Minus][pointIndex] = false;
      borderPointsSubdomainAreValid[6][(int)Mesh::face_t::face1Minus][pointIndex] = false;
    }
  }

  // subdomains 3,7
  pointIndex = 0;
  streamlineIndex--;  // consider center point twice
  for (int i = nBorderPointsX_-1; i < iEndHorizontal; i++, streamlineIndex++, pointIndex++)
  {
    // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
    if (streamlineZPoints[streamlineIndex].size() == nBorderPointsZNew_)
    {

      // loop over bottom half of the streamline points
      for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
      {
        const int subdomainIndex = 3;
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Minus][pointIndex] = true;
      }

      // loop over top half of the streamline points
      for (int zLevelIndex = nBorderPointsZ_-1; zLevelIndex < nBorderPointsZNew_; zLevelIndex++)
      {
        const int subdomainIndex = 7;
        const int zLevelIndexNew = zLevelIndex - (nBorderPointsZ_-1);
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face1Minus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face1Minus][pointIndex] = true;
      }
    }
    else
    {
      borderPointsSubdomainAreValid[3][(int)Mesh::face_t::face1Minus][pointIndex] = false;
      borderPointsSubdomainAreValid[7][(int)Mesh::face_t::face1Minus][pointIndex] = false;
    }
  }

  LOG(DEBUG) << "starting with vertical center line, streamlineIndex = " << streamlineIndex;

  // vertical center line (with corner points and center point)
  // subdomains 0,4, 1,5
  pointIndex = iBeginVertical;
  for (int i = iBeginVertical; i < nBorderPointsX_; i++, streamlineIndex++, pointIndex++)
  {
    // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
    if (streamlineZPoints[streamlineIndex].size() == nBorderPointsZNew_)
    {
      // subdomains 0,4,
      // loop over bottom half of the streamline points
      for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
      {
        const int subdomainIndex = 0;
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Plus][pointIndex] = true;
      }

      // loop over top half of the streamline points
      for (int zLevelIndex = nBorderPointsZ_-1; zLevelIndex < nBorderPointsZNew_; zLevelIndex++)
      {
        const int subdomainIndex = 4;
        const int zLevelIndexNew = zLevelIndex - (nBorderPointsZ_-1);
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Plus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Plus][pointIndex] = true;
      }
    }
    else
    {
      borderPointsSubdomainAreValid[0][(int)Mesh::face_t::face0Plus][pointIndex] = false;
      borderPointsSubdomainAreValid[4][(int)Mesh::face_t::face0Plus][pointIndex] = false;
    }

    // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
    if (streamlineZPoints[streamlineIndex].size() == nBorderPointsZNew_)
    {
      // subdomains 1,5
      // loop over bottom half of the streamline points
      for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
      {
        const int subdomainIndex = 1;
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Minus][pointIndex] = true;
      }

      // loop over top half of the streamline points
      for (int zLevelIndex = nBorderPointsZ_-1; zLevelIndex < nBorderPointsZNew_; zLevelIndex++)
      {
        const int subdomainIndex = 5;
        const int zLevelIndexNew = zLevelIndex - (nBorderPointsZ_-1);
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Minus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Minus][pointIndex] = true;
      }
    }
    else
    {
      borderPointsSubdomainAreValid[1][(int)Mesh::face_t::face0Minus][pointIndex] = false;
      borderPointsSubdomainAreValid[5][(int)Mesh::face_t::face0Minus][pointIndex] = false;
    }
  }

  LOG(DEBUG) << "nBorderPointsZ_: " << nBorderPointsZ_ << ", nBorderPointsX_: " << nBorderPointsX_ << ", iEndVertical: " << iEndVertical;

  // subdomains 2,6, 3,7
  pointIndex = 0;
  streamlineIndex--;  // consider center point twice
  for (int i = nBorderPointsX_-1; i < iEndVertical; i++, streamlineIndex++, pointIndex++)
  {
    // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
    if (streamlineZPoints[streamlineIndex].size() == nBorderPointsZNew_)
    {
      LOG(DEBUG) << "";
      LOG(DEBUG) << "i=" << i << " (iEndVertical=" << iEndVertical << ")";
      LOG(DEBUG) << "streamlineIndex: " << streamlineIndex;
      LOG(DEBUG) << "pointIndex: " << pointIndex << ", nBorderPointsZ_: " << nBorderPointsZ_;

      // subdomains 2,6
      // loop over bottom half of the streamline points
      for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
      {
        const int subdomainIndex = 2;
        LOG(DEBUG) << "i=" << i << " (iEndVertical=" << iEndVertical << "), set borderPointsSubdomain[" << subdomainIndex << "]["
          << (int)Mesh::face_t::face0Plus << "][" << zLevelIndex << "][" <<pointIndex << "] = streamlineZPoints[" << streamlineIndex
          << "][" << zLevelIndex << "] = " << streamlineZPoints[streamlineIndex][zLevelIndex];
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Plus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Plus][pointIndex] = true;
      }

      // loop over top half of the streamline points
      for (int zLevelIndex = nBorderPointsZ_-1; zLevelIndex < nBorderPointsZNew_; zLevelIndex++)
      {
        const int subdomainIndex = 6;
        const int zLevelIndexNew = zLevelIndex - (nBorderPointsZ_-1);
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Plus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Plus][pointIndex] = true;
      }
    }
    else
    {
      borderPointsSubdomainAreValid[2][(int)Mesh::face_t::face0Plus][pointIndex] = false;
      borderPointsSubdomainAreValid[6][(int)Mesh::face_t::face0Plus][pointIndex] = false;
    }

    // check if streamline is valid, i.e. goes from bottom to top, only then use it otherwise approximate streamline from neighbouring streamlines
    if (streamlineZPoints[streamlineIndex].size() == nBorderPointsZNew_)
    {
      // subdomains 3,7
      // loop over bottom half of the streamline points
      for (int zLevelIndex = 0; zLevelIndex < nBorderPointsZ_; zLevelIndex++)
      {
        const int subdomainIndex = 3;
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Minus][zLevelIndex][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Minus][pointIndex] = true;
      }

      // loop over top half of the streamline points
      for (int zLevelIndex = nBorderPointsZ_-1; zLevelIndex < nBorderPointsZNew_; zLevelIndex++)
      {
        const int subdomainIndex = 7;
        const int zLevelIndexNew = zLevelIndex - (nBorderPointsZ_-1);
        borderPointsSubdomain[subdomainIndex][(int)Mesh::face_t::face0Minus][zLevelIndexNew][pointIndex] = streamlineZPoints[streamlineIndex][zLevelIndex];
        borderPointsSubdomainAreValid[subdomainIndex][(int)Mesh::face_t::face0Minus][pointIndex] = true;
      }
    }
    else
    {
      borderPointsSubdomainAreValid[3][(int)Mesh::face_t::face0Minus][pointIndex] = false;
      borderPointsSubdomainAreValid[7][(int)Mesh::face_t::face0Minus][pointIndex] = false;
    }
  }

#ifndef NDEBUG
  std::stringstream s;
  for (int subdomainIndex = 0; subdomainIndex < 8; subdomainIndex++)
  {
    s.str("");
    s << "06_border_points_subdomain_" << subdomainIndex;
    PyObject_CallFunction(functionOutputBorderPoints_, "s i O f", s.str().c_str(), currentRankSubset_->ownRankNo(),
                          PythonUtility::convertToPython<std::array<std::vector<std::vector<Vec3>>,4>>::get(borderPointsSubdomain[subdomainIndex]), 1.0);
    PythonUtility::checkForError();
  }
#endif
}

};  // namespace
