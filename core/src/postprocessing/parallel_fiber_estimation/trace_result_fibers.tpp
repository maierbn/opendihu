#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
bool ParallelFiberEstimation<BasisFunctionType>::
checkTraceFinalFibers()
{
  // determine current level = log2(nRanksPerCoordinateDirection_)
  level_ = 0;
  int nRanksPerCoordinateDirection = nRanksPerCoordinateDirection_[0];
  while (nRanksPerCoordinateDirection >>= 1)
  {
    level_++;
  }

  int nRanksAvailable = this->context_.partitionManager()->nRanksCommWorld();

  // decide if the algorithm should no more refine the subdomains but trace the final fibers
  bool traceFinalFibers = false;
  if (level_ == maxLevel_)
  {
    traceFinalFibers = true;
  }

  if (!traceFinalFibers)
  {
    if (nRanksAvailable < currentRankSubset_->size()*8)
    {
      LOG(WARNING) << "Cannot run algorithm until level " << maxLevel_ << ", currently at level " << level_
        << ", total number of ranks is " << nRanksAvailable << ", number needed for next level would be " << currentRankSubset_->size()*8 << "." << std::endl
        << "Perform final step of algorithm now at level " << level_;
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

      // reorder streamline points such that they go from bottom to top
      if (streamlineDirection < 0 && nRanksZ != 1)
      {
        std::reverse(streamlinePoints.begin(), streamlinePoints.end());
      }

      // sample the streamline at equidistant z levels and store points in fibers vector
      int fibersPointIndex = j * (nFineGridFibers_+1) * nFibersX + i * (nFineGridFibers_+1);
      sampleStreamlineAtEquidistantZPoints(streamlinePoints, seedPoint, bottomZClip, topZClip, fibers[fibersPointIndex], fibersPointIndex*1000);  // last parameter is value for debugging output
    }
  }

  // check which of the key fibers are invalid
  std::vector<std::vector<bool>> keyFiberIsValid(nBorderPointsXNew_, std::vector<bool>(nBorderPointsXNew_, true));  // this only stores for key fibers if they are valid, dimensions are nBorderPointsXNew_ x nBorderPointsXNew_

  int nValid = 0;
  for (int j = 0; j < nBorderPointsXNew_; j++)
  {
    for (int i = 0; i < nBorderPointsXNew_; i++)
    {
      int pointIndex = j * (nFineGridFibers_+1) * nFibersX + i * (nFineGridFibers_+1);

      bool isValid = true;
      if (fibers[pointIndex].size() != nBorderPointsZNew_)
      {
        LOG(DEBUG) << "fiber[" << pointIndex << "] (" << i << "," << j << ") of (" << nBorderPointsXNew_ << "," << nBorderPointsXNew_ << ")"
          << " is not long enough (size: " << fibers[pointIndex].size() << ")";
        isValid = false;
      }
      else if (MathUtility::norm<3>(fibers[pointIndex][0]) < 1e-4)
      {
        LOG(DEBUG) << "fiber[" << pointIndex << "] (" << i << "," << j << ") of (" << nBorderPointsXNew_ << "," << nBorderPointsXNew_ << ")"
          << ", first point is zero";
        isValid = false;
      }
      else if (MathUtility::norm<3>(fibers[pointIndex][1]) < 1e-4)
      {
        LOG(DEBUG) << "fiber[" << pointIndex << "] (" << i << "," << j << ") of (" << nBorderPointsXNew_ << "," << nBorderPointsXNew_ << ")"
          << ", second point is zero";
        isValid = false;
      }
      else if (MathUtility::norm<3>(fibers[pointIndex][nBorderPointsZNew_/2]) < 1e-4)
      {
        LOG(DEBUG) << "fiber[" << pointIndex << "] (" << i << "," << j << ") of (" << nBorderPointsXNew_ << "," << nBorderPointsXNew_ << ")"
          << ", center point is zero";
        isValid = false;
      }
      else
      {
        nValid++;
      }

      if (!isValid)
      {
        keyFiberIsValid[j][i] = false;
      }
    }
  }

  int nInvalid = MathUtility::sqr(nBorderPointsXNew_) - nValid;
  LOG(DEBUG) << "key fibers, number: " << MathUtility::sqr(nBorderPointsXNew_) << ", valid: " << nValid << ", invalid: " << nInvalid;

  // fix the invalid key fibers in the interior by interpolating from the neighbouring fibers
  LOG(DEBUG) << "fixInvalidKeyFibers";
  int nFibersFixed = 0;
  fixInvalidKeyFibers(nFibersX, keyFiberIsValid, fibers, nFibersFixed);

  // send end points of streamlines to next rank that continues the streamline
  exchangeSeedPointsAfterTracingKeyFibers(nRanksZ, rankZNo, nFibersX, streamlineDirectionUpwards, seedPoints, fibers);

  // reduced number of valid/invalid fibers on ranks
  int nValidGlobal = 0;
  int nInvalidGlobal = 0;
  int nFibersFixedGlobal = 0;
  MPI_Reduce(&nValid, &nValidGlobal, 1, MPI_INT, MPI_SUM, 0, currentRankSubset_->mpiCommunicator());
  MPI_Reduce(&nInvalid, &nInvalidGlobal, 1, MPI_INT, MPI_SUM, 0, currentRankSubset_->mpiCommunicator());
  MPI_Reduce(&nFibersFixed, &nFibersFixedGlobal, 1, MPI_INT, MPI_SUM, 0, currentRankSubset_->mpiCommunicator());

  if (currentRankSubset_->ownRankNo() == 0)
  {
    LOG(INFO) << "total number of key fibers, initially valid: " << nValidGlobal << ", initially invalid: " << nInvalidGlobal << ", fixed: " << nFibersFixedGlobal
      << ", finally invalid: " << nInvalidGlobal - nFibersFixedGlobal << ", finally valid: " << nValidGlobal + nFibersFixedGlobal;
  }


#ifndef NDEBUG
#ifdef STL_OUTPUT
  PyObject_CallFunction(functionOutputPoints_, "s i O f", "13_final_seed_points", currentRankSubset_->ownRankNo(),
                        PythonUtility::convertToPython<std::vector<Vec3>>::get(seedPoints), 0.2);
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

  int nFibersGlobal = 0;
  int nFibersNotInterpolatedGlobal = 0;
  MPI_Reduce(&nFibers, &nFibersGlobal, 1, MPI_INT, MPI_SUM, 0, currentRankSubset_->mpiCommunicator());
  MPI_Reduce(&nFibersNotInterpolated, &nFibersNotInterpolatedGlobal, 1, MPI_INT, MPI_SUM, 0, currentRankSubset_->mpiCommunicator());

  if (currentRankSubset_->ownRankNo() == 0)
  {
    LOG(INFO) << "total number of fibers, valid: " << nFibersGlobal-nFibersNotInterpolatedGlobal << ", invalid: " << nFibersNotInterpolatedGlobal;
  }

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
  const int nParameters = 10;   // if I change this, also change the constant header length further down in resampleFibersInFile
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
          fibers[fibersPointIndex].clear();
          fibers[fibersPointIndex].resize(nPointsCurrentFiber, Vec3({0.0,0.0,0.0}));
        }
        if (fibers[fibersPointIndex][int(nPointsCurrentFiber/2)][0] == 0.0 && fibers[fibersPointIndex][int(nPointsCurrentFiber/2)][1] == 0.0 
          && fibers[fibersPointIndex][int(nPointsCurrentFiber/2)][2] == 0.0)
        {
          fibers[fibersPointIndex].clear();
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

    if (nPointsWholeFiber != nNodesPerFiber_)
    {
      resampleFibersInFile(nPointsWholeFiber);
    }
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
    // determine size of file
    struct stat statBuffer;
    stat(resultFilename_.c_str(), &statBuffer);
    int fileSize = statBuffer.st_size;

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
    const long long fiberDataSize = nPointsPerFiber*3*sizeof(double);

    LOG(DEBUG) << "headerLength: " << headerLength << ", nFibers: " << nFibers << ", nPointsPerFiber: " << nPointsPerFiber << ", nFibersX: " << nFibersX;

    if (int((fileSize-(32+headerLength)) / fiberDataSize) != nFibers)
    {
      LOG(ERROR) << "File \"" << resultFilename_ << "\" states to have " << nFibers << " fibers in header, but actually has "
        << int((fileSize-(32+headerLength)) / fiberDataSize) << " fibers!";
    }

    for (int iterationNo = 0; iterationNo < 3 && (nFibersInvalid - nFibersFixed > 0 || iterationNo == 0); iterationNo++)
    {

      // determine which fibers are valid
      std::vector<std::vector<bool>> fiberIsValid(nFibersX, std::vector<bool>(nFibersX, true));

      for (int fiberIndexY = 0; fiberIndexY != nFibersX; fiberIndexY++)
      {
        for (int fiberIndexX = 0; fiberIndexX != nFibersX; fiberIndexX++)
        {
          int fiberIndex = fiberIndexY*nFibersX + fiberIndexX;
          file.seekg(32+headerLength + fiberIndex*fiberDataSize);

          // read all points of fiber
          for (int zPointIndex = 0; zPointIndex < nPointsPerFiber; zPointIndex++)
          {
            Vec3 point;
            MathUtility::readPoint(file, point);

            // if fiber is invalid
            if (point[0] == 0.0 && point[1] == 0.0 && point[2] == 0.0)
            {
              fiberIsValid[fiberIndexY][fiberIndexX] = false;
              nFibersInvalid++;
              LOG(DEBUG) << "fiber (" << fiberIndexX << "," << fiberIndexY << ") / (" << nFibersX << "," << nFibersX << ") is invalid (" << zPointIndex << ". point)";
              break;
            }
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
                assert(previousFiberIndex < nFibers);
                file.seekg(32+headerLength + previousFiberIndex*fiberDataSize + zIndex*3*sizeof(double));

                MathUtility::readPoint(file, point0);

                int nextFiberIndex = fiberIndexY*nFibersX + rightNeighbourIndex;
                assert(nextFiberIndex < nFibers);

                file.seekg(32+headerLength + nextFiberIndex*fiberDataSize + zIndex*3*sizeof(double));

                MathUtility::readPoint(file, point1);

                // compute the interpolated value
                Vec3 interpolatedPoint = alpha0*point0 + alpha1*point1;

                // write the interpolated value back
                int interpolatedFiberIndex = fiberIndexY*nFibersX + fiberIndexX;
                assert(interpolatedFiberIndex < nFibers);

                file.seekp(32+headerLength + interpolatedFiberIndex*fiberDataSize + zIndex*3*sizeof(double));

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
                assert(previousFiberIndex < nFibers);

                file.seekg(32+headerLength + previousFiberIndex*fiberDataSize + zIndex*3*sizeof(double));

                MathUtility::readPoint(file, point0);

                int nextFiberIndex = backNeighbourIndex*nFibersX + fiberIndexX;
                assert(nextFiberIndex < nFibers);

                file.seekg(32+headerLength + nextFiberIndex*fiberDataSize + zIndex*3*sizeof(double));

                MathUtility::readPoint(file, point1);

                // compute the interpolated value
                Vec3 interpolatedPoint = alpha0*point0 + alpha1*point1;

                // write the interpolated value back
                int interpolatedFiberIndex = fiberIndexY*nFibersX + fiberIndexX;
                assert(interpolatedFiberIndex < nFibers);

                file.seekp(32+headerLength + interpolatedFiberIndex*fiberDataSize + zIndex*3*sizeof(double));

                MathUtility::writePoint(file, interpolatedPoint);
              }
              nFibersFixed++;
            }
            else
            {
              // check if first point of invalid fiber is set
              int fiberIndex = fiberIndexY*nFibersX + fiberIndexX;
              file.seekg(32+headerLength + fiberIndex*fiberDataSize);

              // read first point of fiber
              Vec3 firstPoint;
              MathUtility::readPoint(file, firstPoint);

              // if fiber is valid
              if (firstPoint[0] != 0.0 || firstPoint[1] != 0.0 || firstPoint[2] != 0.0)
              {
                int validIndex0 = -1;
                int validIndex1 = -1;
                if (rightNeighbourIndex != -1)
                {
                  validIndex0 = fiberIndexY*nFibersX + rightNeighbourIndex;
                }
                else if (leftNeighbourIndex != -1)
                {
                  validIndex0 = fiberIndexY*nFibersX + leftNeighbourIndex;
                }

                if (backNeighbourIndex != -1)
                {
                  validIndex1 = backNeighbourIndex*nFibersX + fiberIndexX;
                }
                else if(frontNeighbourIndex != -1)
                {
                  validIndex1 = frontNeighbourIndex*nFibersX + fiberIndexX;
                }

                if (validIndex0 != -1 && validIndex1 != -1)
                {

                  // read first points of fibers

                  // read first point of fiber
                  std::array<Vec3,2> neighbouringFiberFirstPoint;

                  file.seekg(32+headerLength + validIndex0*fiberDataSize);
                  MathUtility::readPoint(file, neighbouringFiberFirstPoint[0]);

                  file.seekg(32+headerLength + validIndex1*fiberDataSize);
                  MathUtility::readPoint(file, neighbouringFiberFirstPoint[1]);

                  // assert that fiber points are on same z position
                  assert(neighbouringFiberFirstPoint[0][2] - neighbouringFiberFirstPoint[1][2] < 1e-12);
                  if (firstPoint[2] - neighbouringFiberFirstPoint[1][2] > 1e-12)
                  {
                    LOG(WARNING) << "Could not fixe fiber (" << fiberIndexX << "," << fiberIndexY << ") / (" << nFibersX << "," << nFibersX << "),"
                     << " first point: " << firstPoint << ", neighbouring points: " << neighbouringFiberFirstPoint[0] << "," << neighbouringFiberFirstPoint[1];

                    int zIndex = nPointsPerFiber - 1;
                    file.seekg(32+headerLength + validIndex0*fiberDataSize + zIndex*3*sizeof(double));
                    MathUtility::readPoint(file, neighbouringFiberFirstPoint[0]);

                    file.seekg(32+headerLength + validIndex1*fiberDataSize + zIndex*3*sizeof(double));
                    MathUtility::readPoint(file, neighbouringFiberFirstPoint[1]);
                  }

                  if (firstPoint[2] - neighbouringFiberFirstPoint[1][2] > 1e-12)
                  {
                    LOG(WARNING) << "Could not fixe fiber (" << fiberIndexX << "," << fiberIndexY << ") / (" << nFibersX << "," << nFibersX << "),"
                     << " first point: " << firstPoint << ", neighbouring points: " << neighbouringFiberFirstPoint[0] << "," << neighbouringFiberFirstPoint[1];

                    int zIndex = int(nPointsPerFiber/2);
                    file.seekg(32+headerLength + validIndex0*fiberDataSize + zIndex*3*sizeof(double));
                    MathUtility::readPoint(file, neighbouringFiberFirstPoint[0]);

                    file.seekg(32+headerLength + validIndex1*fiberDataSize + zIndex*3*sizeof(double));
                    MathUtility::readPoint(file, neighbouringFiberFirstPoint[1]);
                  }

                  if (firstPoint[2] - neighbouringFiberFirstPoint[1][2] > 1e-12)
                  {
                    LOG(ERROR) << "Could not fixed fiber (" << fiberIndexX << "," << fiberIndexY << ") / (" << nFibersX << "," << nFibersX << "),"
                      << " first point: " << firstPoint << ", neighbouring points: " << neighbouringFiberFirstPoint[0] << "," << neighbouringFiberFirstPoint[1];
                  }
                  else
                  {
                    //assert(firstPoint[2] - neighbouringFiberFirstPoint[1][2] < 1e-12);

                    // the two neighbouring valid fibers are the indices validIndex0 and validIndex1
                    // consider the triangle between the first valid point is of the to be interpolated fiber and the two neighbouring fibers
                    double angle = std::atan2(firstPoint[1]-neighbouringFiberFirstPoint[0][1], firstPoint[0]-neighbouringFiberFirstPoint[0][0]);
                    double relativeLength = MathUtility::distance<3>(firstPoint, neighbouringFiberFirstPoint[0])
                      / MathUtility::distance<3>(neighbouringFiberFirstPoint[1], neighbouringFiberFirstPoint[0]);

                    // interpolate points
                    // loop over all points of the fiber
                    for (int zIndex = 0; zIndex != nPointsPerFiber; zIndex++)
                    {
                      // get the two valid fibers
                      Vec3 point0, point1;
                      assert(validIndex0 < nFibers);
                      file.seekg(32+headerLength + validIndex0*fiberDataSize + zIndex*3*sizeof(double));

                      MathUtility::readPoint(file, point0);

                      assert(validIndex1 < nFibers);
                      file.seekg(32+headerLength + validIndex1*fiberDataSize + zIndex*3*sizeof(double));

                      MathUtility::readPoint(file, point1);

                      double distance = MathUtility::distance<3>(point0, point1);

                      // compute the interpolated point
                      Vec3 interpolatedPoint = point0 + Vec3({cos(angle)*relativeLength*distance, sin(angle)*relativeLength*distance, 0.0});

                      // write the interpolated value back
                      int interpolatedFiberIndex = fiberIndexY*nFibersX + fiberIndexX;
                      assert(interpolatedFiberIndex < nFibers);

                      file.seekp(32+headerLength + interpolatedFiberIndex*fiberDataSize + zIndex*3*sizeof(double));

                      MathUtility::writePoint(file, interpolatedPoint);
                    }
                    nFibersFixed++;
                  }
                }
              }
            }
          }
        }
      }

      LOG(DEBUG) << "nFibersInvalid: " << nFibersInvalid << ", nFibersFixed: " << nFibersFixed << ", difference: " << nFibersInvalid - nFibersFixed;
      LOG(INFO) << "Iteration " << iterationNo << ", fixed " << nFibersFixed << " of " << nFibersInvalid << " invalid fibers, total: " << nFibers;
    }
  }

  file.close();
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
resampleFibersInFile(int nPointsPerFiber)
{
  // create a new file with all the fibers from the old file but resampled such that they have nNodesPerFiber_ nodes
  LOG(INFO) << "resample fibers in file";

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

  const int headerLength = sizeof(int32_t)*10;

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

  double oldZIncrement = double(topZClip_ - bottomZClip_) / (nPointsPerFiber - 1);

  int nFibersX = (nBorderPointsXNew_-1) * nFineGridFibers_ + nBorderPointsXNew_;
  int nFibers = MathUtility::sqr(nFibersX);
  const long long fiberDataSize = nPointsPerFiber*3*sizeof(double);

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

      //LOG(DEBUG) << "clip: " << bottomZClip_ << "," << topZClip_ << " zIndex: " << zIndex << "/" << nNodesPerFiber_ << ", currentZ: " << currentZ;

      // compute z indices of the point in the old fiber between which the new point will be
      int oldZIndexPrevious = int((currentZ-bottomZClip_) / oldZIncrement);
      int oldZIndexNext = oldZIndexPrevious + 1;

      // load previous point, only if it was not already loaded
      if (oldZIndexPrevious != oldZIndexPrevious0)
      {
        fileOld.seekg(32+headerLength + fiberIndex*fiberDataSize + oldZIndexPrevious*3*sizeof(double));
        MathUtility::readPoint(fileOld, previousPoint);
        oldZIndexPrevious0 = oldZIndexPrevious;
      }

      // load next point
      if (oldZIndexNext < nPointsPerFiber)
      {
        // only if it was not already loaded
        if (oldZIndexNext != oldZIndexNext0)
        {
          fileOld.seekg(32+headerLength + fiberIndex*fiberDataSize + oldZIndexNext*3*sizeof(double));
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
        LOG(DEBUG) << "f" << fiberIndex << " z" << zIndex << "(" << currentZ << ") indices " << oldZIndexPrevious << "," << oldZIndexNext
          << ", points " << previousPoint << nextPoint << ", alpha: " << alpha << ", newPoint: " << newPoint;
      }

      // write point to file
      MathUtility::writePoint(fileNew, newPoint);
    }
  }

  fileOld.close();
  fileNew.close();
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
fixInvalidKeyFibers(int nFibersX, std::vector<std::vector<bool>> &fiberIsValid, std::vector<std::vector<Vec3>> &fibers, int &nFibersFixed)
{
  LOG(DEBUG) << "fixInvalidFibers";

  // fibers[fiberIndex][zLevelIndex]
  nFibersFixed = 0;

  // loop over invalid fibers and fix them from neighbouring fibers
  for (int fiberIndexY = 0; fiberIndexY != nBorderPointsXNew_; fiberIndexY++)
  {
    for (int fiberIndexX = 0; fiberIndexX != nBorderPointsXNew_; fiberIndexX++)
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
        for (int neighbouringFiberIndexX = fiberIndexX+1; neighbouringFiberIndexX < nBorderPointsXNew_; neighbouringFiberIndexX++)
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
        for (int neighbouringFiberIndexY = fiberIndexY+1; neighbouringFiberIndexY < nBorderPointsXNew_; neighbouringFiberIndexY++)
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
          for (int zIndex = 0; zIndex != nBorderPointsZNew_; zIndex++)
          {
            // get the left and right valid fibers
            int previousFiberIndex = (fiberIndexY*nFibersX + leftNeighbourIndex) * (nFineGridFibers_+1);
            int nextFiberIndex = (fiberIndexY*nFibersX + rightNeighbourIndex) * (nFineGridFibers_+1);

            Vec3 point0 = fibers[previousFiberIndex][zIndex];
            Vec3 point1 = fibers[nextFiberIndex][zIndex];

            // compute the interpolated value
            Vec3 interpolatedPoint = alpha0*point0 + alpha1*point1;

            // write the interpolated value back
            int interpolatedFiberIndex = (fiberIndexY*nFibersX + fiberIndexX) * (nFineGridFibers_+1);
            fibers[interpolatedFiberIndex][zIndex] = interpolatedPoint;
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
          for (int zIndex = 0; zIndex != nBorderPointsZNew_; zIndex++)
          {
            // get the left and right valid fibers
            int previousFiberIndex = (frontNeighbourIndex*nFibersX + fiberIndexX) * (nFineGridFibers_+1);
            int nextFiberIndex = (backNeighbourIndex*nFibersX + fiberIndexX) * (nFineGridFibers_+1);

            Vec3 point0 = fibers[previousFiberIndex][zIndex];
            Vec3 point1 = fibers[nextFiberIndex][zIndex];

            // compute the interpolated value
            Vec3 interpolatedPoint = alpha0*point0 + alpha1*point1;

            // write the interpolated value back
            int interpolatedFiberIndex = (fiberIndexY*nFibersX + fiberIndexX) * (nFineGridFibers_+1);
            fibers[interpolatedFiberIndex][zIndex] = interpolatedPoint;
          }
          nFibersFixed++;
        }
      }
    }
  }

  LOG(DEBUG) << "n key fibers fixed: " << nFibersFixed;
}

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
interpolateFineFibersFromFile()
{
  // open existing file to read

  struct stat statBuffer;
  int rc = stat(resultFilename_.c_str(), &statBuffer);
  if (rc != 0)
  {
    LOG(FATAL) << "Could not open file \"" << resultFilename_ << "\".";
  }
  int fileOldSize = statBuffer.st_size;

  std::ifstream fileOld(resultFilename_.c_str(), std::ios::in | std::ios::binary);
  assert (fileOld.is_open());

  std::stringstream newFilename;

  // check if old file name is like "7x7fibers"
  bool filenameHasXFormat = false;
  bool xFound = false;
  int suffixPos = 0;
  for(int i = 0; i < resultFilename_.size(); i++)
  {
    if (!isdigit(resultFilename_[i]))
    {
      if (xFound)
      {
        suffixPos = i;
        filenameHasXFormat = true;
        break;
      }
      if (resultFilename_[i] == 'x')
      {
        xFound = true;
      }
      else
      {
        break;
      }
    }
  }

  // parse header
  // skip first part of header
  fileOld.seekg(32);
  union int32
  {
    char c[4];
    int32_t i;
  }
  bufferHeaderLength, bufferNFibers, bufferNPointsPerFiber;

  // get length of header
  fileOld.read(bufferHeaderLength.c, 4);
  int headerLength = bufferHeaderLength.i;
  assert(headerLength == sizeof(int32_t)*10);

  // get number of fibers
  fileOld.read(bufferNFibers.c, 4);
  int nFibersOld = bufferNFibers.i;

  // get number of points per fiber
  fileOld.read(bufferNPointsPerFiber.c, 4);
  int nPointsPerFiber = bufferNPointsPerFiber.i;

  int nFibersOldX = int(std::round(std::sqrt(nFibersOld)));

  LOG(DEBUG) << "read header: nFibers: " << nFibersOld << ", nFibersX: " << nFibersOldX << ", nPointsPerFiber: " << nPointsPerFiber;

  // reset to beginning of header
  fileOld.seekg(0);

  // copy header
  std::vector<char> headerBuffer(32+headerLength);
  fileOld.read(headerBuffer.data(), 32+headerLength);

  int nFibersNewX = (nFibersOldX-1) * nFineGridFibers_ + nFibersOldX;
  int nFibersNew = MathUtility::sqr(nFibersNewX);

  if (filenameHasXFormat)
  {
    newFilename << nFibersNewX << "x" << resultFilename_.substr(suffixPos);
  }
  else
  {
    newFilename << resultFilename_ << ".fine";
  }

  LOG(DEBUG) << "write to file " << newFilename.str();
  std::ofstream fileNew(newFilename.str().c_str(), std::ios::out | std::ios::binary);
  if (!fileNew.is_open())
  {
    LOG(FATAL) << "Could not write to file \"" << newFilename.str() << "\".";
  }

  // write header to new file
  fileNew.write(headerBuffer.data(), 32+headerLength);

  LOG(INFO) << "\"" << resultFilename_ << "\" -> \"" << newFilename.str() << "\"";
  LOG(INFO) << "interpolate from " <<  nFibersOldX << " x " << nFibersOldX << " = " << nFibersOld
    << " to " << nFibersNewX << " x " << nFibersNewX << " = " << nFibersNew << ", (nFineGridFibers: " << nFineGridFibers_ << ")";

  const long long fiberDataSize = nPointsPerFiber*3*sizeof(double);
  long long nBytes = (32+headerLength) + (long long)(nFibersNew)*fiberDataSize;
  long long nGibiBytes = nBytes / (long long)(1024) / (long long)(1024) / (long long)(1024);
  LOG(INFO) << "estimated size: " << nGibiBytes << " GiB";

  if (nGibiBytes >= 1)
  {
    LOG(INFO) << "Press any key to continue . . .";
    std::cin.get();
  }

  if (int((fileOldSize-(32+headerLength)) / fiberDataSize) != nFibersOld)
  {
    LOG(ERROR) << "File \"" << resultFilename_ << "\" states to have " << nFibersOld << " fibers in header, but actually has "
      << int((fileOldSize-(32+headerLength)) / fiberDataSize) << " fibers!";
  }

  union
  {
    int32_t parameter;
    char c[sizeof(int32_t)];
  };

  // write new number of fibers
  // nFibersTotal
  fileNew.seekp(32+4);
  parameter = nFibersNew;
  fileNew.write(c, 4);

  // nBorderPointsXNew
  fileNew.seekp(32+3*4);
  parameter = nFibersOldX;
  fileNew.write(c, 4);

  // nFineGridFibers
  fileNew.seekp(32+5*4);
  parameter = nFineGridFibers_;
  fileNew.write(c, 4);

  // date
  fileNew.seekp(32+9*4);
  parameter = time(NULL);
  fileNew.write(c, 4);

  fileNew.seekp(32+headerLength);

  // define 4 buffers that can each hold data for one fiber
  std::array<std::vector<char>,4> buffer;

  // initialize buffers
  for (int i = 0; i < 4; i++)
  {
    buffer[i].resize(fiberDataSize);
  }

  // loop over fibers in new file
  for (int fiberIndexY = 0; fiberIndexY != nFibersNewX; fiberIndexY++)
  {
    for (int fiberIndexX = 0; fiberIndexX != nFibersNewX; fiberIndexX++)
    {
      int oldFiberIndexX = (int)(fiberIndexX / (nFineGridFibers_+1));
      int oldFiberIndexY = (int)(fiberIndexY / (nFineGridFibers_+1));

      VLOG(1) << "fiber " << fiberIndexX << "," << fiberIndexY << ", oldFiberIndex: " << oldFiberIndexX << "," << oldFiberIndexY;

      // if fiber is key fiber
      if (fiberIndexX % (nFineGridFibers_+1) == 0 && fiberIndexY % (nFineGridFibers_+1) == 0)
      {
        // copy fiber data from old file
        int oldFiberIndex = oldFiberIndexY * nFibersOldX + oldFiberIndexX;

        // assert that file is long enough
        if (32+headerLength + (oldFiberIndex+1)*fiberDataSize > fileOldSize)
        {
          LOG(ERROR) << "Read past end of file, oldFiberIndex: " << oldFiberIndex;
        }
        assert(32+headerLength + (oldFiberIndex+1)*fiberDataSize <= fileOldSize);

        fileOld.seekg(32+headerLength + oldFiberIndex*fiberDataSize);
        fileOld.read(buffer[0].data(), fiberDataSize);
        fileNew.write(buffer[0].data(), fiberDataSize);
      }
      else
      {
        // if fiber is no key fiber, read 4 neighbouring fibers and interpolate
        // determine indices
        int oldFiberIndex[4] = {
          oldFiberIndexY * nFibersOldX + oldFiberIndexX,
          oldFiberIndexY * nFibersOldX + (oldFiberIndexX+1),
          (oldFiberIndexY+1) * nFibersOldX + oldFiberIndexX,
          (oldFiberIndexY+1) * nFibersOldX + (oldFiberIndexX+1)
        };

        std::array<bool,4> edgeFiberIsValid({true,true,true,true});
        for (int i = 0; i < 4; i++)
        {
          std::size_t fileEndPos = 32+headerLength + (oldFiberIndex[i]+1)*fiberDataSize;
          if (fileEndPos > fileOldSize)
          {
            edgeFiberIsValid[i] = false;
            continue;
          }

          fileOld.seekg(32+headerLength + oldFiberIndex[i]*fiberDataSize);
          fileOld.read(buffer[i].data(), fiberDataSize);
        }

        // compute interpolation factors
        double alpha0 = double(fiberIndexX - oldFiberIndexX*(nFineGridFibers_+1)) / (nFineGridFibers_+1);
        double alpha1 = double(fiberIndexY - oldFiberIndexY*(nFineGridFibers_+1)) / (nFineGridFibers_+1);

        // 2 3
        // 0 1
        if (!edgeFiberIsValid[2] || !edgeFiberIsValid[3])
        {
          alpha1 = 0.0;
        }

        VLOG(1) << "(" << fiberIndexX << "," << fiberIndexY << ") alpha: " << alpha0 << "," << alpha1;

        // loop over points of fiber
        for (int zPointIndex = 0; zPointIndex < nPointsPerFiber; zPointIndex++)
        {
          // get edge points
          std::array<Vec3,4> point;
          for (int i = 0; i < 4; i++)
          {
            for (int j = 0; j < 3; j++)
            {
              double *memoryLocation = reinterpret_cast<double *>(&buffer[i][(zPointIndex*3 + j)*sizeof(double)]);
              point[i][j] = *memoryLocation;
            }
          }

          Vec3 resultPoint
            = (1.-alpha0) * (1.-alpha1) * point[0]
            + alpha0     * (1.-alpha1) * point[1]
            + (1.-alpha0) * alpha1     * point[2]
            + alpha0     * alpha1     * point[3];

          // write result point
          fileNew.write((char *)resultPoint.data(), 3*sizeof(double));
        }
      }
    }
  }

  fileNew.close();
}

} // namespace
