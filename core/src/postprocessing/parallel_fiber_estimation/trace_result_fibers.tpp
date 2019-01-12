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
        << ", total number of ranks is " << nRanksAvailable << ", number needed for next level would be " << currentRankSubset_->size()*8 << "." << std::endl
        << "Perform final step of algorithm now at level " << level;
      traceFinalFibers = true;
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
  pointIndexStride = nFibersX * (nFineGridFibers_+1);
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
        VLOG(1) << "set value at left border, pointIndexStride: " << pointIndexStride << ", fibers[" << fibersPointIndex << "][z"
          << zLevelIndex << "] =" << fibers[fibersPointIndex][zLevelIndex];
        VLOG(1) << "set value at left border, pointIndexStride: " << pointIndexStride << ", fibers[" << fibersPointIndex << "][z"
          << zLevelIndex+(nBorderPointsZ_-1) << "] =" << fibers[fibersPointIndex][zLevelIndex+(nBorderPointsZ_-1)];
      }
    }
  }

  // right (0+)
  subdomainIndices = std::vector<int>{1, 3};
  face = Mesh::face_t::face0Plus;
  pointIndexStart = nFibersX-1;
  pointIndexStride = nFibersX * (nFineGridFibers_+1);
  fibersPointIndex = pointIndexStart;
  for (int i = 0; i < subdomainIndices.size(); i++)
  {
    const int subdomainIndex = subdomainIndices[i];
    for (int pointIndex = 0; pointIndex != nBorderPointsX_-(i==1? 0: 1); pointIndex++, fibersPointIndex += pointIndexStride)
    {
      for (int zLevelIndex = 0; zLevelIndex != nBorderPointsZ_; zLevelIndex++)
      {
        fibers[fibersPointIndex][zLevelIndex] = borderPointsSubdomain[subdomainIndex][face][zLevelIndex][pointIndex];
        fibers[fibersPointIndex][zLevelIndex+(nBorderPointsZ_-1)] = borderPointsSubdomain[subdomainIndex+4][face][zLevelIndex][pointIndex];
        VLOG(1) << "set value at right border, pointIndexStride: " << pointIndexStride << ", fibers[" << fibersPointIndex << "][z"
          << zLevelIndex << "] =" << fibers[fibersPointIndex][zLevelIndex];
        VLOG(1) << "set value at right border, pointIndexStride: " << pointIndexStride << ", fibers[" << fibersPointIndex << "][z"
          << zLevelIndex+(nBorderPointsZ_-1) << "] =" << fibers[fibersPointIndex][zLevelIndex+(nBorderPointsZ_-1)];
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

      // for serial execution, the seed point is at the center of the streamline
      if (nRanksZ == 1)
      {
        // trace streamlines forwards
        std::vector<Vec3> forwardPoints;
        this->traceStreamline(seedPoint, 1.0, forwardPoints);

        // trace streamline backwards
        std::vector<Vec3> backwardPoints;
        this->traceStreamline(seedPoint, -1.0, backwardPoints);

        // copy collected points to result vector, note avoiding this additional copy-step is not really possible, since it would require a push_front which is only efficient with lists, but we need a vector here
        streamlinePoints.insert(streamlinePoints.begin(), backwardPoints.rbegin(), backwardPoints.rend());
        streamlinePoints.insert(streamlinePoints.end(), seedPoint);
        streamlinePoints.insert(streamlinePoints.end(), forwardPoints.begin(), forwardPoints.end());

      }
      else
      {
        // with parallel execution the seed point is at one end of the fiber
        streamlinePoints.push_back(seedPoint);

        this->traceStreamline(seedPoint, streamlineDirection, streamlinePoints);
      }

      // if everything was cleared, add seed point
      if (streamlinePoints.empty())
        streamlinePoints.push_back(seedPoint);

      // save end points of streamlines such that neighbouring rank can continue
      streamlineEndPoints[j*nBorderPointsXNew_+i][0] = streamlinePoints.back();

      // reorder streamline points such that they go from bottom to top
      if (streamlineDirection < 0 && nRanksZ != 1)
      {
        std::reverse(streamlinePoints.begin(), streamlinePoints.end());
      }

      // sample the streamline at equidistant z levels and store points in fibers vector
      int fibersPointIndex = j * (nFineGridFibers_+1) * nFibersX + i * (nFineGridFibers_+1);
      sampleStreamlineAtEquidistantZPoints(streamlinePoints, seedPoint, bottomZClip, topZClip, fibers[fibersPointIndex], fibersPointIndex*1000);
    }
  }

  // check which fibers are invalid
  int nValid = 0;
  for (int j = 0; j < nBorderPointsXNew_; j++)
  {
    for (int i = 0; i < nBorderPointsXNew_; i++)
    {
      int pointIndex = j * (nFineGridFibers_+1) * nFibersX + i * (nFineGridFibers_+1);

      if (fibers[pointIndex].size() != nBorderPointsZNew_)
      {
        LOG(DEBUG) << "fiber[" << pointIndex << "] (" << i << "," << j << ") is not long enough (size: " << fibers[pointIndex].size() << ")";
      }
      else if (MathUtility::norm<3>(fibers[pointIndex][0]) < 1e-4)
      {
        LOG(DEBUG) << "fiber[" << pointIndex << "] (" << i << "," << j << "), first point is zero";
      }
      else if (MathUtility::norm<3>(fibers[pointIndex][1]) < 1e-4)
      {
        LOG(DEBUG) << "fiber[" << pointIndex << "] (" << i << "," << j << "), second point is zero";
      }
      else if (MathUtility::norm<3>(fibers[pointIndex][nBorderPointsZNew_/2]) < 1e-4)
      {
        LOG(DEBUG) << "fiber[" << pointIndex << "] (" << i << "," << j << "), center point is zero";
      }
      else
      {
        nValid++;
      }
    }
  }

  LOG(DEBUG) << "key fibers, number: " << MathUtility::sqr(nBorderPointsXNew_) << ", valid: " << nValid << ", invalid: " << MathUtility::sqr(nBorderPointsXNew_) - nValid;

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

      // determine if key fibers are all valid
      bool keyFibersInvalid = false;
      if (fibers[keyFiberPointIndex0].size() != nBorderPointsZNew_ || fibers[keyFiberPointIndex1].size() != nBorderPointsZNew_ || fibers[keyFiberPointIndex2].size() != nBorderPointsZNew_ || fibers[keyFiberPointIndex3].size() != nBorderPointsZNew_)
      {
        LOG(DEBUG) << "fibers at (" << i << "," << j << ") are incomplete, sizes: " << fibers[keyFiberPointIndex0].size()
          << "," << fibers[keyFiberPointIndex1].size()
          << "," << fibers[keyFiberPointIndex2].size()
          << "," << fibers[keyFiberPointIndex3].size();
        nFibersNotInterpolated += MathUtility::sqr(nFineGridFibers_+1);
        keyFibersInvalid = true;
      }

      if (!keyFibersInvalid)
      {
        if (MathUtility::norm<3>(fibers[keyFiberPointIndex0][nBorderPointsZNew_/2]) < 1e-4 || MathUtility::norm<3>(fibers[keyFiberPointIndex1][nBorderPointsZNew_/2]) < 1e-4
          || MathUtility::norm<3>(fibers[keyFiberPointIndex2][nBorderPointsZNew_/2]) < 1e-4 || MathUtility::norm<3>(fibers[keyFiberPointIndex3][nBorderPointsZNew_/2]) < 1e-4)
        {
          LOG(DEBUG) << "fibers at (" << i << "," << j << ") are invalid: " << std::endl
            << fibers[keyFiberPointIndex0] << std::endl
            << fibers[keyFiberPointIndex1] << std::endl
            << fibers[keyFiberPointIndex2] << std::endl
            << fibers[keyFiberPointIndex3];
          keyFibersInvalid = true;
        }
      }

      if (!keyFibersInvalid)
      {
        VLOG(1) << "interpolate fine fibers at (" << i << "," << j << "), index " << keyFiberPointIndex0;
      }

      for (int fineGridJ = 0; fineGridJ < nFineGridFibers_+1; fineGridJ++)
      {
        for (int fineGridI = 0; fineGridI < nFineGridFibers_+1; fineGridI++)
        {
          // do not consider key fiber
          if (fineGridJ == 0 && fineGridI == 0)
            continue;

          // get index of fiber in fibers vector
          int fibersPointIndex = keyFiberPointIndex0 + fineGridJ*nFibersX + fineGridI;

          // if key fibers are invalid, set interpolated fibers to 0
          if (keyFibersInvalid)
          {
            for (int zLevelIndex = 0; zLevelIndex != nBorderPointsZNew_; zLevelIndex++)
            {
              fibers[fibersPointIndex][zLevelIndex] = Vec3({0.0,0.0,0.0});
            }
          }
          else
          {
            // interpolate fibers using key fibers
            Vec2 xi({fineGridI / (nFineGridFibers_ + 1.0), fineGridJ / (nFineGridFibers_ + 1.0)});

            VLOG(1) << "  fine fiber (" << fineGridI << "," << fineGridJ << ") index " << fibersPointIndex << ", xi: " << xi;

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
  }
#if 0
  PyObject_CallFunction(functionOutputStreamlines_, "s i O f", "14_final", currentRankSubset_->ownRankNo(),
                        PythonUtility::convertToPython<std::vector<std::vector<Vec3>>>::get(fibers), 0.1);
  PythonUtility::checkForError();
#endif
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
  const int nParameters = 10;
  if (ownRankNo == 0)
  {
    std::string writeBuffer("opendihu binary fibers file     ");

    union
    {
      int32_t parameters[nParameters];
      char c[nParameters*sizeof(int32_t)];
    };

    parameters[0] = nParameters*sizeof(int32_t);
    parameters[1] = nFibersTotal;
    parameters[2] = nPointsWholeFiber;
    parameters[3] = nBorderPointsXNew_;
    parameters[4] = nBorderPointsZNew_;
    parameters[5] = nFineGridFibers_;
    parameters[6] = currentRankSubset_->size();
    parameters[7] = nRanksZ;
    parameters[8] = nFibers;
    parameters[9] = time(NULL);

    LOG(DEBUG) << "nFibersTotal: " << nFibersTotal << ", nPointsWholeFiber: " << nPointsWholeFiber << ", nFibersPerRank: " << nFibers;

    writeBuffer += std::string(c, nParameters*sizeof(int32_t));

    MPI_Status status;
    MPIUtility::handleReturnValue(MPI_File_write(fileHandle, writeBuffer.c_str(), writeBuffer.length(), MPI_BYTE, &status), "MPI_File_write", &status);

    headerOffset = writeBuffer.size();
    LOG(DEBUG) << "headerOffset: " << headerOffset;
  }

  // set hardcoded header offset
  assert(headerOffset == 32+nParameters*sizeof(int32_t) || ownRankNo != 0);
  headerOffset = 32+nParameters*sizeof(int32_t);

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

        // if fibers are invalid, fill with 0.0
        if (fibers[fibersPointIndex].size() < nPointsCurrentFiber)
        {
          fibers[fibersPointIndex].resize(nPointsCurrentFiber, Vec3({0.0,0.0,0.0}));
        }

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

      VLOG(1) << "write fiber (" << i << "," << j << "), global (" << fiberIndex0 << "," << fiberIndex1 << "), global index "
        << fiberIndex << ", nValues: " << nValues << ", nBytes: " << nBytes
        << ", pointOffset: " << pointOffset << " (" << nPointsWholeFiber << " points per fiber), headerOffset: " << headerOffset << ", offset: " << offset;

      // write to file
      MPI_Status status;
      MPIUtility::handleReturnValue(MPI_File_write_at(fileHandle, offset, writeBuffer, nBytes, MPI_BYTE, &status), "MPI_File_write_at", &status);
    }
  }

  MPIUtility::handleReturnValue(MPI_File_close(&fileHandle), "MPI_File_close");

  LOG(DEBUG) << "Fibers written to file \"" << resultFilename_ << "\".";

  if (currentRankSubset_->ownRankNo() == 0)
  {
    fixInvalidFibersInFile();
  }
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
fixInvalidFibersInFile()
{
  // open the file again and interpolate all missing fibers

  int nPointsPerFiber = 0;
  int nFibers = 0;
  int headerLength = 0;
  std::fstream file(resultFilename_.c_str(), std::ios::out | std::ios::in | std::ios::binary);
  if (!file.is_open())
  {
    LOG(FATAL) << "Could not open file \"" << resultFilename_ << "\".";
  }
  else
  {
    // skip first part of header
    file.seekg(32);
    union int32
    {
      char c[4];
      int32_t i;
    }
    bufferHeaderLength, bufferNFibers, bufferNPointsPerFiber;

    // get length of header
    file.read(bufferHeaderLength.c, 4);
    headerLength = bufferHeaderLength.i;

    // get number of fibers
    file.read(bufferNFibers.c, 4);
    nFibers = bufferNFibers.i;

    // get number of points per fiber
    file.read(bufferNPointsPerFiber.c, 4);
    nPointsPerFiber = bufferNPointsPerFiber.i;

    // skip rest of header
    file.seekg(32+headerLength);

    int nFibersX = int(std::round(std::sqrt(nFibers)));
    int nFibersInvalid = 0;
    int nFibersFixed = 0;

    LOG(DEBUG) << "headerLength: " << headerLength << ", nFibers: " << nFibers << ", nPointsPerFiber: " << nPointsPerFiber << ", nFibersX: " << nFibersX;

    // determine which fibers are valid
    std::vector<std::vector<bool>> fiberIsValid(nFibersX, std::vector<bool>(nFibers, true));

    for (int fiberIndexY = 0; fiberIndexY != nFibersX; fiberIndexY++)
    {
      for (int fiberIndexX = 0; fiberIndexX != nFibersX; fiberIndexX++)
      {
        int fiberIndex = fiberIndexY*nFibersX + fiberIndexX;
        file.seekg(32+headerLength + fiberIndex*nPointsPerFiber*3*sizeof(double));

        // read first value of fiber
        Vec3 firstPoint;
        MathUtility::readPoint(file, firstPoint);

        // if fiber is invalid
        if (firstPoint[0] == 0.0 && firstPoint[1] == 0.0 && firstPoint[2] == 0.0)
        {
          fiberIsValid[fiberIndexY][fiberIndexX] = false;
          nFibersInvalid++;
        }
      }
    }

    VLOG(2) << "fiberIsValid: " << fiberIsValid;

    // loop over invalid fibers and fix them from neighbouring fibers
    for (int fiberIndexY = 0; fiberIndexY != nFibersX; fiberIndexY++)
    {
      for (int fiberIndexX = 0; fiberIndexX != nFibersX; fiberIndexX++)
      {
        if (!fiberIsValid[fiberIndexY][fiberIndexX])
        {
          // find neighbouring valid fibers
          std::vector<std::pair<int,int>> neighbouringFibers;

          // check X direction
          // find previous neighbour
          int leftNeighbourIndex = -1;
          for (int neighbouringFiberIndexX = fiberIndexX-1; neighbouringFiberIndexX >= 0; neighbouringFiberIndexX--)
          {
            if (fiberIsValid[fiberIndexY][neighbouringFiberIndexX])
            {
              leftNeighbourIndex = neighbouringFiberIndexX;
              break;
            }
          }

          // find next neighoubr
          int rightNeighbourIndex = -1;
          for (int neighbouringFiberIndexX = fiberIndexX+1; neighbouringFiberIndexX < nFibersX; neighbouringFiberIndexX++)
          {
            if (fiberIsValid[fiberIndexY][neighbouringFiberIndexX])
            {
              rightNeighbourIndex = neighbouringFiberIndexX;
              break;
            }
          }

          // check Y direction
          // find previous neighbour
          int frontNeighbourIndex = -1;
          for (int neighbouringFiberIndexY = fiberIndexY-1; neighbouringFiberIndexY >= 0; neighbouringFiberIndexY--)
          {
            if (fiberIsValid[neighbouringFiberIndexY][fiberIndexX])
            {
              frontNeighbourIndex = neighbouringFiberIndexY;
              break;
            }
          }

          // find next neighoubr
          int backNeighbourIndex = -1;
          for (int neighbouringFiberIndexY = fiberIndexY+1; neighbouringFiberIndexY < nFibersX; neighbouringFiberIndexY++)
          {
            if (fiberIsValid[neighbouringFiberIndexY][fiberIndexX])
            {
              backNeighbourIndex = neighbouringFiberIndexY;
              break;
            }
          }

          VLOG(2) << "fiber " << fiberIndexX << "," << fiberIndexY << ": neighbours x("
            << leftNeighbourIndex << "," << rightNeighbourIndex << ") y(" << frontNeighbourIndex << "," << backNeighbourIndex << ")";

          if (leftNeighbourIndex != -1 && rightNeighbourIndex != -1)
          {
            // interpolate fiber from the neighbours in x direction

            double alpha0 = double(rightNeighbourIndex - fiberIndexX) / (rightNeighbourIndex - leftNeighbourIndex);
            double alpha1 = double(fiberIndexX - leftNeighbourIndex) / (rightNeighbourIndex - leftNeighbourIndex);

            VLOG(2) << "take left-right, alphas: " << alpha0 << ", " << alpha1;

            // loop over all points of the fiber
            for (int zIndex = 0; zIndex != nPointsPerFiber; zIndex++)
            {
              // get the left and right valid fibers
              Vec3 point0, point1;
              int previousFiberIndex = fiberIndexY*nFibersX + leftNeighbourIndex;
              file.seekg(32+headerLength + previousFiberIndex*nPointsPerFiber*3*sizeof(double) + zIndex*3*sizeof(double));

              MathUtility::readPoint(file, point0);

              int nextFiberIndex = fiberIndexY*nFibersX + rightNeighbourIndex;
              file.seekg(32+headerLength + nextFiberIndex*nPointsPerFiber*3*sizeof(double) + zIndex*3*sizeof(double));

              MathUtility::readPoint(file, point1);

              // compute the interpolated value
              Vec3 interpolatedPoint = alpha0*point0 + alpha1*point1;

              // write the interpolated value back
              int interpolatedFiberIndex = fiberIndexY*nFibersX + fiberIndexX;
              file.seekp(32+headerLength + interpolatedFiberIndex*nPointsPerFiber*3*sizeof(double) + zIndex*3*sizeof(double));

              MathUtility::writePoint(file, interpolatedPoint);
            }
            nFibersFixed++;
          }
          else if (frontNeighbourIndex != -1 && backNeighbourIndex != -1)
          {
            // interpolate fiber from the neighbours in x direction

            double alpha0 = double(backNeighbourIndex - fiberIndexX) / (backNeighbourIndex - frontNeighbourIndex);
            double alpha1 = double(fiberIndexX - frontNeighbourIndex) / (backNeighbourIndex - frontNeighbourIndex);

            VLOG(2) << "take front-back, alphas: " << alpha0 << ", " << alpha1;

            // loop over all points of the fiber
            for (int zIndex = 0; zIndex != nPointsPerFiber; zIndex++)
            {
              // get the left and right valid fibers
              Vec3 point0, point1;
              int previousFiberIndex = frontNeighbourIndex*nFibersX + fiberIndexX;
              file.seekg(32+headerLength + previousFiberIndex*nPointsPerFiber*3*sizeof(double) + zIndex*3*sizeof(double));

              MathUtility::readPoint(file, point0);

              int nextFiberIndex = backNeighbourIndex*nFibersX + fiberIndexX;
              file.seekg(32+headerLength + nextFiberIndex*nPointsPerFiber*3*sizeof(double) + zIndex*3*sizeof(double));

              MathUtility::readPoint(file, point1);

              // compute the interpolated value
              Vec3 interpolatedPoint = alpha0*point0 + alpha1*point1;

              // write the interpolated value back
              int interpolatedFiberIndex = fiberIndexY*nFibersX + fiberIndexX;
              file.seekp(32+headerLength + interpolatedFiberIndex*nPointsPerFiber*3*sizeof(double) + zIndex*3*sizeof(double));

              MathUtility::writePoint(file, interpolatedPoint);
            }
            nFibersFixed++;
          }
        }
      }
    }

    LOG(DEBUG) << "nFibersInvalid: " << nFibersInvalid << ", nFibersFixed: " << nFibersFixed << ", difference: " << nFibersInvalid - nFibersFixed;
  }

  file.close();

  if (nPointsPerFiber != nNodesPerFiber_)
  {
    // create a new file with all the fibers from the old file but resampled such that they have nNodesPerFiber_ nodes

    // rename existing file
    std::string newFilename = resultFilename_ + std::string("_");
    std::stringstream moveCommand;
    moveCommand << "mv " << resultFilename_ << " " << newFilename;
    int ret = std::system(moveCommand.str().c_str());
    ret++;
    std::this_thread::sleep_for(std::chrono::milliseconds(100));

    // open existing file to read
    std::ifstream fileOld(newFilename.c_str(), std::ios::in | std::ios::binary);
    assert (fileOld.is_open());

    LOG(DEBUG) << "write to file " << resultFilename_;
    std::ofstream fileNew(resultFilename_.c_str(), std::ios::out | std::ios::binary);
    assert (fileNew.is_open());

    // copy header
    std::vector<char> headerBuffer(32+headerLength);
    fileOld.read(headerBuffer.data(), 32+headerLength);
    fileNew.write(headerBuffer.data(), 32+headerLength);

    // write new number of fibers
    union
    {
      int32_t parameter;
      char c[sizeof(int32_t)];
    };
    parameter = nNodesPerFiber_;

    fileNew.seekp(32+8);
    fileNew.write(c, 4);

    fileNew.seekp(32+headerLength);

    double oldZIncrement = bottomZClip_ + double(topZClip_ - bottomZClip_) / (nPointsPerFiber - 1);

    // loop over fibers in old file
    for (int fiberIndex = 0; fiberIndex != nFibers; fiberIndex++)
    {
      Vec3 previousPoint;
      Vec3 nextPoint;

      int oldZIndexPrevious0 = -1;    /// z index of the point that is currently loadad in previousPoint
      int oldZIndexNext0 = -1;        /// z index of the point that is currently loadad in nextPoint

      // loop over nodes of the new fiber and write them to the new file
      for (int zIndex = 0; zIndex < nNodesPerFiber_; zIndex++)
      {
        // compute the z value of the current point in the new fiber
        double currentZ = bottomZClip_ + zIndex * double(topZClip_ - bottomZClip_) / (nNodesPerFiber_ - 1);

        // compute z indices of the point in the old fiber between which the new point will be
        int oldZIndexPrevious = int((currentZ-bottomZClip_) / oldZIncrement);
        int oldZIndexNext = oldZIndexPrevious + 1;

        // load previous point, only if it was not already loaded
        if (oldZIndexPrevious != oldZIndexPrevious0)
        {
          fileOld.seekg(32+headerLength + fiberIndex*nPointsPerFiber*3*sizeof(double) + oldZIndexPrevious*3*sizeof(double));
          MathUtility::readPoint(fileOld, previousPoint);
          oldZIndexPrevious0 = oldZIndexPrevious;
        }

        // load next point
        if (oldZIndexNext < nPointsPerFiber)
        {
          // only if it was not already loaded
          if (oldZIndexNext != oldZIndexNext0)
          {
            fileOld.seekg(32+headerLength + fiberIndex*nPointsPerFiber*3*sizeof(double) + oldZIndexNext*3*sizeof(double));
            MathUtility::readPoint(fileOld, nextPoint);
            oldZIndexNext0 = oldZIndexNext;
          }
        }
        else
        {
          nextPoint = previousPoint;
        }

        // compute new point by interpolation between previousPoint and nextPoint
        double alpha = (currentZ - (bottomZClip_ + oldZIndexPrevious*oldZIncrement)) / oldZIncrement;
        Vec3 newPoint = (1.-alpha) * previousPoint + alpha * nextPoint;

        if (fiberIndex < 10 || fiberIndex > nFibers-10)
        {
          VLOG(2) << "f" << fiberIndex << " z" << zIndex << "(" << currentZ << ") indices " << oldZIndexPrevious << "," << oldZIndexNext
            << ", points " << previousPoint << nextPoint << ", alpha: " << alpha << ", newPoint: " << newPoint;
        }

        // write point to file
        MathUtility::writePoint(fileNew, newPoint);
      }
    }

    fileOld.close();
    fileNew.close();
  }
}

} // namespace
