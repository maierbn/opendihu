#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
bool ParallelFiberEstimation<BasisFunctionType>::
checkTraceFinalFibers(int &level)
{
  // determine current level = log2(nRanksPerCoordinateDirection_)
  level = 0;
  int nRanksPerCoordinateDirection = nRanksPerCoordinateDirection_[0];
  while (nRanksPerCoordinateDirection >>= 1)
  {
    level++;
  }

  int nRanksAvailable = this->context_.partitionManager()->nRanksCommWorld();

  // decide if the algorithm should no more refine the subdomains but trace the final fibers
  bool traceFinalFibers = false;
  if (level == maxLevel_)
  {
    traceFinalFibers = true;
  }

  if (!traceFinalFibers)
  {
    if (nRanksAvailable < currentRankSubset_->size()*8)
    {
      LOG(WARNING) << "Cannot run algorithm until level " << maxLevel_ << ", currently at level " << level
        << ", total number of ranks is " << nRanksAvailable << ", number needed for next level would be " << currentRankSubset_->size()*8 << ".";
      //traceFinalFibers = true;
    }
  }
  return traceFinalFibers;
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
traceResultFibers(double streamlineDirection, int seedPointsZIndex, const std::vector<Vec3> &nodePositions,
                  std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &borderPointsSubdomain)
{
  LOG(DEBUG) << "traceResultFibers, nBorderPointsXNew_: " << nBorderPointsXNew_ << ", nFineGridFibers: " << nFineGridFibers_;

  // grid, x = key fibers which are traced, o = interpolated fibers for nFineGridFibers_ = 1
  //  _ _ _ _
  // |x o x o x|
  // |o o o o o|
  // |x o x o x|
  // |o o o o o|
  // |x_o_x_o_x|

  // total number of fibers per coordinate direction, including the ones on the left and right borders
  int nFibersX = (nBorderPointsXNew_-1) * nFineGridFibers_ + nBorderPointsXNew_;
  //int nFibersX = nBorderPointsXNew_*(nFineGridFibers_+1) - nFineGridFibers_;

  int nFibers = MathUtility::sqr(nFibersX);

  // allocate memory
  std::vector<std::vector<Vec3>> fibers(nFibers);     // [fiberIndex][zLevelIndex]
  for (int fiberIndex = 0; fiberIndex != nFibers; fiberIndex++)
  {
    fibers[fiberIndex].resize(nBorderPointsZNew_,Vec3{0,0,0});
  }

  // fill already existing fibers at borders

  //   ^ --(1+)-> ^   ^ --(1+)-> ^
  //   0-   [2]   0+  0-   [3]   0+
  //   | --(1-)-> |   | --(1-)-> |
  //
  //   ^ --(1+)-> ^   ^ --(1+)-> ^
  // ^ 0-   [0]   0+  0-   [1]   0+
  // | | --(1-)-> |   | --(1-)-> |
  // +-->

  // bottom (1-)
  std::vector<int> subdomainIndices = {0, 1};
  int face = Mesh::face_t::face1Minus;
  int pointIndexStart = 0;
  int pointIndexStride = nFineGridFibers_ + 1;
  int fibersPointIndex = pointIndexStart;
  for (int i = 0; i < subdomainIndices.size(); i++)
  {
    const int subdomainIndex = subdomainIndices[i];
    for (int pointIndex = 0; pointIndex != nBorderPointsX_-1; pointIndex++, fibersPointIndex += pointIndexStride)
    {
      for (int zLevelIndex = 0; zLevelIndex != nBorderPointsZ_; zLevelIndex++)
      {
        fibers[fibersPointIndex][zLevelIndex] = borderPointsSubdomain[subdomainIndex][face][zLevelIndex][pointIndex];
        fibers[fibersPointIndex][zLevelIndex+(nBorderPointsZ_-1)] = borderPointsSubdomain[subdomainIndex+4][face][zLevelIndex][pointIndex];
      }
    }
  }

  // top (1+)
  subdomainIndices = std::vector<int>{2, 3};
  face = Mesh::face_t::face1Plus;
  pointIndexStart = nFibersX * (nFibersX-1);
  fibersPointIndex = pointIndexStart;
  for (int i = 0; i < subdomainIndices.size(); i++)
  {
    const int subdomainIndex = subdomainIndices[i];
    for (int pointIndex = 0; pointIndex != nBorderPointsX_-1; pointIndex++, fibersPointIndex += pointIndexStride)
    {
      for (int zLevelIndex = 0; zLevelIndex != nBorderPointsZ_; zLevelIndex++)
      {
        fibers[fibersPointIndex][zLevelIndex] = borderPointsSubdomain[subdomainIndex][face][zLevelIndex][pointIndex];
        fibers[fibersPointIndex][zLevelIndex+(nBorderPointsZ_-1)] = borderPointsSubdomain[subdomainIndex+4][face][zLevelIndex][pointIndex];
      }
    }
  }

  // left (0-)
  subdomainIndices = std::vector<int>{0, 2};
  face = Mesh::face_t::face0Minus;
  pointIndexStart = 0;
  pointIndexStride = nFibersX;
  fibersPointIndex = pointIndexStart;
  for (int i = 0; i < subdomainIndices.size(); i++)
  {
    const int subdomainIndex = subdomainIndices[i];
    for (int pointIndex = 0; pointIndex != nBorderPointsX_-1; pointIndex++, fibersPointIndex += pointIndexStride)
    {
      for (int zLevelIndex = 0; zLevelIndex != nBorderPointsZ_; zLevelIndex++)
      {
        fibers[fibersPointIndex][zLevelIndex] = borderPointsSubdomain[subdomainIndex][face][zLevelIndex][pointIndex];
        fibers[fibersPointIndex][zLevelIndex+(nBorderPointsZ_-1)] = borderPointsSubdomain[subdomainIndex+4][face][zLevelIndex][pointIndex];
      }
    }
  }

  // right (0+)
  subdomainIndices = std::vector<int>{1, 3};
  face = Mesh::face_t::face0Plus;
  pointIndexStart = nFibersX-1;
  pointIndexStride = nFibersX;
  fibersPointIndex = pointIndexStart;
  for (int i = 0; i < subdomainIndices.size(); i++)
  {
    const int subdomainIndex = subdomainIndices[i];
    for (int pointIndex = 0; pointIndex != nBorderPointsX_-1; pointIndex++, fibersPointIndex += pointIndexStride)
    {
      for (int zLevelIndex = 0; zLevelIndex != nBorderPointsZ_; zLevelIndex++)
      {
        fibers[fibersPointIndex][zLevelIndex] = borderPointsSubdomain[subdomainIndex][face][zLevelIndex][pointIndex];
        fibers[fibersPointIndex][zLevelIndex+(nBorderPointsZ_-1)] = borderPointsSubdomain[subdomainIndex+4][face][zLevelIndex][pointIndex];
      }
    }
  }

#ifndef NDEBUG
  PyObject_CallFunction(functionOutputStreamlines_, "s i O f", "12_final_border", currentRankSubset_->ownRankNo(),
                        PythonUtility::convertToPython<std::vector<std::vector<Vec3>>>::get(fibers), 0.1);
  PythonUtility::checkForError();
#endif

  // determine own seed points
  // set seed points for interior fibers and trace fibers
  std::vector<Vec3> seedPoints(nBorderPointsXNew_*nBorderPointsXNew_, Vec3({0.0,0.0,0.0}));
  for (int j = 1; j < nBorderPointsXNew_-1; j++)
  {
    for (int i = 1; i < nBorderPointsXNew_-1; i++)
    {
      // determine seed point
      Vec3 seedPoint = nodePositions[seedPointsZIndex*nBorderPointsXNew_*nBorderPointsXNew_ + j*nBorderPointsXNew_ + i];
      seedPoints[j*nBorderPointsXNew_+i] = seedPoint;
    }
  }

  // communicate seed Points
  int rankZNo = meshPartition_->ownRankPartitioningIndex(2);
  int nRanksZ = meshPartition_->nRanks(2);
  bool streamlineDirectionUpwards = streamlineDirection>0;

  // The overall picture is that global streamlines begin at the center (at rankZNo/2).
  // Rank int(nRanksZ/2) send the initial seed points to the rank below (int(nRanksZ/2)-1)
  // Then every rank traces its streamlines and sends the end points as new seed points to the next rank (lower or upper neighbour, depending on streamlineDirection)

  // determine if previously set seedPoints are used or if they are received from neighbouring rank, receive seed points or send them to lower neighbour, if own rank is int(nRanksZ/2)
  exchangeSeedPointsBeforeTracing(nRanksZ, rankZNo, streamlineDirectionUpwards, seedPoints);

  // determine z range of current subdomain
  double bottomZClip = 0;
  double topZClip = 0;
  computeBottomTopZClip(bottomZClip, topZClip);

  // set seed points for interior fibers and trace fibers
  // define and initialize vector that contains the end points of the traced streamlines.
  // These serve as new seed points at the neighbouring rank which continues the streamlines.
  std::vector<std::vector<Vec3>> streamlineEndPoints(nBorderPointsXNew_*nBorderPointsXNew_);
  for (int i = 0; i < nBorderPointsXNew_*nBorderPointsXNew_; i++)
    streamlineEndPoints[i].resize(1,Vec3({0.0,0.0,0.0}));

  for (int j = 1; j < nBorderPointsXNew_-1; j++)
  {
    for (int i = 1; i < nBorderPointsXNew_-1; i++)
    {
      // determine seed point
      Vec3 seedPoint = seedPoints[j*nBorderPointsXNew_+i];

      // trace streamline
      std::vector<Vec3> streamlinePoints;
      streamlinePoints.push_back(seedPoint);

      this->traceStreamline(seedPoint, streamlineDirection, streamlinePoints);

      // if everything was cleared, add seed point
      if (streamlinePoints.empty())
        streamlinePoints.push_back(seedPoint);

      // save end points of streamlines such that neighbouring rank can continue
      streamlineEndPoints[j*nBorderPointsXNew_+i][0] = streamlinePoints.back();

      // reorder streamline points such that they go from bottom to top
      if (streamlineDirection < 0)
      {
        std::reverse(streamlinePoints.begin(), streamlinePoints.end());
      }

      // sample the streamline at equidistant z levels and store points in fibers vector
      int fibersPointIndex = j * (nFineGridFibers_+1) * nFibersX + i * (nFineGridFibers_+1);
      sampleStreamlineAtEquidistantZPoints(streamlinePoints, seedPoint, bottomZClip, topZClip, fibers[fibersPointIndex], fibersPointIndex*1000);
    }
  }

  // send end points of streamlines to next rank that continues the streamline
  exchangeSeedPointsAfterTracing(nRanksZ, rankZNo, streamlineDirectionUpwards, seedPoints, streamlineEndPoints);

#ifndef NDEBUG
#ifdef STL_OUTPUT
  PyObject_CallFunction(functionOutputPoints_, "s i O f", "13_final_seed_points", currentRankSubset_->ownRankNo(),
                        PythonUtility::convertToPython<std::vector<Vec3>>::get(seedPoints), 0.2);
  PythonUtility::checkForError();
#endif
#endif


#ifndef NDEBUG
#ifdef STL_OUTPUT
  PyObject_CallFunction(functionOutputStreamlines_, "s i O f", "13_streamline_end_points", currentRankSubset_->ownRankNo(),
                        PythonUtility::convertToPython<std::vector<std::vector<Vec3>>>::get(streamlineEndPoints), 0.2);
  PythonUtility::checkForError();
#endif
#endif

#ifndef NDEBUG
  PyObject_CallFunction(functionOutputStreamlines_, "s i O f", "13_final_interior", currentRankSubset_->ownRankNo(),
                        PythonUtility::convertToPython<std::vector<std::vector<Vec3>>>::get(fibers), 0.1);
  PythonUtility::checkForError();
#endif

  // interpolate fibers in between
  LOG(DEBUG) << "interpolate fibers in between key fibers, nFineGridFibers_: " << nFineGridFibers_;
  int nFibersNotInterpolated = 0;

  // grid, x = key fibers which are traced, o = interpolated fibers for nFineGridFibers_ = 1
  //  _ _ _ _
  // |x o x o x|
  // |o o o o o|
  // |x o x o x|
  // |o o o o o|
  // |x_o_x_o_x|
  for (int j = 0; j < nBorderPointsXNew_-1; j++)
  {
    for (int i = 0; i < nBorderPointsXNew_-1; i++)
    {
      // |x   x
      // |o o
      // |x_o_x
      int keyFiberPointIndex0 = j * (nFineGridFibers_+1) * nFibersX + i * (nFineGridFibers_+1);
      int keyFiberPointIndex1 = j * (nFineGridFibers_+1) * nFibersX + (i+1) * (nFineGridFibers_+1);
      int keyFiberPointIndex2 = (j+1) * (nFineGridFibers_+1) * nFibersX + i * (nFineGridFibers_+1);
      int keyFiberPointIndex3 = (j+1) * (nFineGridFibers_+1) * nFibersX + (i+1) * (nFineGridFibers_+1);

      if (fibers[keyFiberPointIndex0].size() != nBorderPointsZNew_ || fibers[keyFiberPointIndex1].size() != nBorderPointsZNew_ || fibers[keyFiberPointIndex2].size() != nBorderPointsZNew_ || fibers[keyFiberPointIndex3].size() != nBorderPointsZNew_)
      {
        LOG(DEBUG) << "fibers at (" << i << "," << j << ") are incomplete, sizes: " << fibers[keyFiberPointIndex0].size()
          << "," << fibers[keyFiberPointIndex1].size()
          << "," << fibers[keyFiberPointIndex2].size()
          << "," << fibers[keyFiberPointIndex3].size();
        nFibersNotInterpolated += MathUtility::sqr(nFineGridFibers_+1);
        continue;
      }
      else
      {
        LOG(DEBUG) << "interpolate fine fibers at (" << i << "," << j << "), index " << keyFiberPointIndex0;
      }

      /*std::array<Vec3,4> quadrilateralPoints;
      quadrilateralPoints[0] = fibers[keyFiberPointIndex0][seedPointsZIndex];
      quadrilateralPoints[1] = fibers[keyFiberPointIndex1][seedPointsZIndex];
      quadrilateralPoints[2] = fibers[keyFiberPointIndex2][seedPointsZIndex];
      quadrilateralPoints[3] = fibers[keyFiberPointIndex3][seedPointsZIndex];*/

      for (int fineGridJ = 0; fineGridJ < nFineGridFibers_+1; fineGridJ++)
      {
        for (int fineGridI = 0; fineGridI < nFineGridFibers_+1; fineGridI++)
        {
          // do not consider key fiber
          if (fineGridJ == 0 && fineGridI == 0)
            continue;

          // get index of fiber in fibers vector
          int fibersPointIndex = keyFiberPointIndex0 + fineGridJ*nFibersX + fineGridI;


          Vec2 xi({fineGridI / (nFineGridFibers_ + 1.0), fineGridJ / (nFineGridFibers_ + 1.0)});

          //Vec3 seedPoint = fibers[fibersPointIndex][seedPointsZIndex];
          // get baricentric coordinates of point seedPoint is quadrilateral
          //MathUtility::quadrilateralGetPointCoordinates(quadrilateralPoints, seedPoint, xi);

          LOG(DEBUG) << "  fine fiber (" << fineGridI << "," << fineGridJ << ") index " << fibersPointIndex << ", xi: " << xi;

          for (int zLevelIndex = 0; zLevelIndex != nBorderPointsZNew_; zLevelIndex++)
          {
            // barycentric interpolation
            fibers[fibersPointIndex][zLevelIndex] =
                (1.0 - xi[0]) * (1.0 - xi[1]) * fibers[keyFiberPointIndex0][zLevelIndex]
              + xi[0]         * (1.0 - xi[1]) * fibers[keyFiberPointIndex1][zLevelIndex]
              + (1.0 - xi[0]) * xi[1]         * fibers[keyFiberPointIndex2][zLevelIndex]
              + xi[0]         * xi[1]         * fibers[keyFiberPointIndex3][zLevelIndex];
          }
        }
      }
    }
  }

  PyObject_CallFunction(functionOutputStreamlines_, "s i O f", "14_final", currentRankSubset_->ownRankNo(),
                        PythonUtility::convertToPython<std::vector<std::vector<Vec3>>>::get(fibers), 0.1);
  PythonUtility::checkForError();

  LOG(DEBUG) << "invalid fibers: " << nFibersNotInterpolated << ", valid fibers: " << nFibers - nFibersNotInterpolated;

  // write fiber data to file, the (x+,x-,y+,y-) border fibers of the global borders are not written to the file because they
  // were not traced but estimated from the border mesh which may be of bad quality
  int ownRankNo = currentRankSubset_->ownRankNo();

  LOG(DEBUG) << "open file \"" << resultFilename_ << "\".";
  // open file
  MPI_File fileHandle;
  MPIUtility::handleReturnValue(MPI_File_open(currentRankSubset_->mpiCommunicator(), resultFilename_.c_str(),
                                              //MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_UNIQUE_OPEN,
                                              MPI_MODE_WRONLY | MPI_MODE_CREATE,
                                              MPI_INFO_NULL, &fileHandle), "MPI_File_open");

  // write file header
  int nPointsWholeFiber = meshPartition_->nRanks(2) * (nBorderPointsZNew_-1) + 1;
  int nFibersRow0 = meshPartition_->nRanks(0) * (nFibersX-1) - 1;
  int nFibersTotal = MathUtility::sqr(nFibersRow0);

  int headerOffset = 0;
  if (ownRankNo == 0)
  {
    std::string writeBuffer("opendihu binary fibers file     ");

    const int nParameters = 8;
    union
    {
      int32_t parameters[nParameters];
      char c[nParameters*sizeof(int32_t)];
    };

    parameters[0] = nFibersTotal;
    parameters[1] = nPointsWholeFiber;
    parameters[2] = nBorderPointsXNew_;
    parameters[3] = nBorderPointsZNew_;
    parameters[4] = nFineGridFibers_;
    parameters[5] = currentRankSubset_->size();
    parameters[6] = nRanksZ;
    parameters[7] = nFibers;

    LOG(DEBUG) << "nFibersTotal: " << nFibersTotal << ", nPointsWholeFiber: " << nPointsWholeFiber << ", nFibersPerRank: " << nFibers;

    writeBuffer += std::string(c, nParameters*sizeof(int32_t));

    MPI_Status status;
    MPIUtility::handleReturnValue(MPI_File_write(fileHandle, writeBuffer.c_str(), writeBuffer.length(), MPI_BYTE, &status), "MPI_File_write", &status);

    headerOffset = writeBuffer.size();
    LOG(DEBUG) << "headerOffset: " << headerOffset;
  }

  // set hardcoded header offset
  assert(headerOffset == 64 || ownRankNo != 0);
  headerOffset = 64;

  // write fibers
  LOG(DEBUG) << "write fibers, nFibersX: " << nFibersX << ", nRanks: (" << meshPartition_->nRanks(0) << "," << meshPartition_->nRanks(1) << "," << meshPartition_->nRanks(2) << ")";
  for (int j = 0; j < nFibersX-1; j++)
  {
    // do not consider first fiber which of a border fiber
    if (j == 0 && meshPartition_->ownRankPartitioningIndex(1) == 0)
      continue;

    for (int i = 0; i < nFibersX-1; i++)
    {
      // only consider last fiber in right column
      if (i == 0 && meshPartition_->ownRankPartitioningIndex(0) == 0)
        continue;

      // convert streamline data to bytes
      int nPointsCurrentFiber = nBorderPointsZNew_ - 1;
      if (meshPartition_->ownRankPartitioningIndex(2) == meshPartition_->nRanks(2)-1)
      {
        nPointsCurrentFiber = nBorderPointsZNew_;
      }

      int nValues = nPointsCurrentFiber*3;
      int nBytes = nValues*sizeof(double);

      std::vector<double> values(nValues);

      // fill values vector
      for (int zLevelIndex = 0; zLevelIndex < nPointsCurrentFiber; zLevelIndex++)
      {
        int fibersPointIndex = j * nFibersX + i;

        for (int k = 0; k < 3; k++)
        {
          values[zLevelIndex*3 + k] = fibers[fibersPointIndex][zLevelIndex][k];
        }
      }
      char *writeBuffer = reinterpret_cast<char *>(values.data());

      // compute offset in file
      int fiberIndex0 = meshPartition_->ownRankPartitioningIndex(0) * (nFibersX-1) + i - 1;
      int fiberIndex1 = meshPartition_->ownRankPartitioningIndex(1) * (nFibersX-1) + j - 1;

      int fiberIndex = fiberIndex1 * nFibersRow0 + fiberIndex0;
      int pointOffset = fiberIndex * nPointsWholeFiber + meshPartition_->ownRankPartitioningIndex(2) * (nBorderPointsZNew_-1);

      int offset = headerOffset + pointOffset * 3 * sizeof(double);

      LOG(DEBUG) << "write fiber (" << i << "," << j << "), global (" << fiberIndex0 << "," << fiberIndex1 << "), global index "
        << fiberIndex << ", nValues: " << nValues << ", nBytes: " << nBytes
        << ", pointOffset: " << pointOffset << " (" << nPointsWholeFiber << " points per fiber), headerOffset: " << headerOffset << ", offset: " << offset;

      // write to file
      MPI_Status status;
      MPIUtility::handleReturnValue(MPI_File_write_at(fileHandle, offset, writeBuffer, nBytes, MPI_BYTE, &status), "MPI_File_write_at", &status);
    }
  }

  MPIUtility::handleReturnValue(MPI_File_close(&fileHandle), "MPI_File_close");

  LOG(DEBUG) << "Fibers written to file \"" << resultFilename_ << "\".";
}

};  // namespace
