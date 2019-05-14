#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
determineLevel()
{
  // determine current level = log2(nRanksPerCoordinateDirection_)
  level_ = 0;
  int nRanksPerCoordinateDirection = nRanksPerCoordinateDirection_[0];
  while (nRanksPerCoordinateDirection >>= 1)
  {
    level_++;
  }
}

template<typename BasisFunctionType>
bool ParallelFiberEstimation<BasisFunctionType>::
checkTraceFinalFibers()
{
  determineLevel();

  int nRanksAvailable = this->context_.nRanksCommWorld();

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
                  const std::array<std::array<std::vector<std::vector<Vec3>>,4>,8> &borderPointsSubdomain, bool finalFile)
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
  PyObject_CallFunction(functionOutputStreamlines_, "s i i O f", "12_final_border", currentRankSubset_->ownRankNo(), level_,
                        PythonUtility::convertToPython<std::vector<std::vector<Vec3>>>::get(fibers), 0.1);
  PythonUtility::checkForError();
#endif

  // determine own seed points
  std::vector<Vec3> seedPoints(nBorderPointsXNew_*nBorderPointsXNew_, Vec3({0.0,0.0,0.0}));

  // set seed points for interior fibers and trace fibers
  for (int j = 1; j < nBorderPointsXNew_-1; j++)
  {
    for (int i = 1; i < nBorderPointsXNew_-1; i++)
    {
      // determine seed point
      Vec3 seedPoint = nodePositions[seedPointsZIndex*nBorderPointsXNew_*nBorderPointsXNew_ + j*nBorderPointsXNew_ + i];
      seedPoints[j*nBorderPointsXNew_+i] = seedPoint;
    }
  }

  LOG(DEBUG) << " set seed points at border";

  // set seed points for border
  // bottom (1-)
  fibersPointIndex = 0;
  pointIndexStride = nFineGridFibers_+1;
  for (int i = 0; i < nBorderPointsXNew_; i++, fibersPointIndex += pointIndexStride)
  {
    seedPoints[0*nBorderPointsXNew_+i] = fibers[fibersPointIndex][seedPointsZIndex];
  }

  // top (1+)
  fibersPointIndex = nFibersX * (nFibersX-1);
  pointIndexStride = nFineGridFibers_+1;
  for (int i = 0; i < nBorderPointsXNew_; i++, fibersPointIndex += pointIndexStride)
  {
    seedPoints[(nBorderPointsXNew_-1)*nBorderPointsXNew_+i] = fibers[fibersPointIndex][seedPointsZIndex];
  }

  // left (0-)
  fibersPointIndex = 0;
  pointIndexStride = nFibersX * (nFineGridFibers_+1);
  for (int i = 0; i < nBorderPointsXNew_; i++, fibersPointIndex += pointIndexStride)
  {
    seedPoints[i*nBorderPointsXNew_+0] = fibers[fibersPointIndex][seedPointsZIndex];
  }

  // right (0+)
  fibersPointIndex = nFibersX-1;
  pointIndexStride = nFibersX * (nFineGridFibers_+1);
  for (int i = 0; i < nBorderPointsXNew_; i++, fibersPointIndex += pointIndexStride)
  {
    seedPoints[i*nBorderPointsXNew_+(nBorderPointsXNew_-1)] = fibers[fibersPointIndex][seedPointsZIndex];
  }

  LOG(DEBUG) << "determine " << nBorderPointsXNew_ << "x" << nBorderPointsXNew_ << " seed points: " << seedPoints.size();

  // communicate seed Points
  int rankZNo = meshPartition_->ownRankPartitioningIndex(2);
  int nRanksZ = meshPartition_->nRanks(2);
  bool streamlineDirectionUpwards = streamlineDirection>0;

  // The overall picture is that global streamlines begin at the center (at rankZNo/2).
  // Rank int(nRanksZ/2) send the initial seed points to the rank below (int(nRanksZ/2)-1)
  // Then every rank traces its streamlines and sends the end points as new seed points to the next rank (lower or upper neighbour, depending on streamlineDirection)

  // determine if previously set seedPoints are used or if they are received from neighbouring rank, receive seed points or send them to lower neighbour, if own rank is int(nRanksZ/2)
  exchangeSeedPointsBeforeTracing(nRanksZ, rankZNo, streamlineDirectionUpwards, seedPoints);

  LOG(DEBUG) << "number of seed points after exchange with neighbours: " << seedPoints.size();

  // determine z range of current subdomain
  double bottomZClip = 0;
  double topZClip = 0;
  computeBottomTopZClip(bottomZClip, topZClip);

  std::vector<std::vector<Vec3>> rawSampledStreamlinesForDebugging;

  // get seed points for interior fibers and trace fibers, the border fibers are already traced
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
      sampleStreamlineAtEquidistantZPoints(streamlinePoints, seedPoint, bottomZClip, topZClip, fibers[fibersPointIndex], fibersPointIndex*1000, rawSampledStreamlinesForDebugging);  // last parameter is value for debugging output

      LOG(DEBUG) << " (" << i << "," << j << ") seedPoint " << seedPoint << ", sampled streamline (" << i << "," << j << ") at equidistant z points";
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
  exchangeSeedPointsAfterTracingKeyFibers(nRanksZ, rankZNo, streamlineDirectionUpwards, nFibersX, seedPoints, fibers);

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
  PyObject_CallFunction(functionOutputPoints_, "s i i O f", "13_final_seed_points", currentRankSubset_->ownRankNo(), level_,
                        PythonUtility::convertToPython<std::vector<Vec3>>::get(seedPoints), 0.2);
  PythonUtility::checkForError();
#endif
#endif

#ifndef NDEBUG
  PyObject_CallFunction(functionOutputStreamlines_, "s i i O f", "13_final_interior", currentRankSubset_->ownRankNo(), level_,
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
  PyObject_CallFunction(functionOutputStreamlines_, "s i i O f", "14_final", currentRankSubset_->ownRankNo(), level_,
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

  int nPointsWholeFiber = meshPartition_->nRanks(2) * (nBorderPointsZNew_-1) + 1;
  int nFibersRow0 = meshPartition_->nRanks(0) * (nFibersX-1) - 1;
  int nFibersTotal = MathUtility::sqr(nFibersRow0);

  std::stringstream filenameStr;
  std::string filename;
  if (finalFile)
  {
    filenameStr << resultFilename_;
  }
  else
  {
    filenameStr << "out/level_" << level_ << "/result_0x0_level_" << level_ << ".bin";
  }

  filename = filenameStr.str();
  // add the number of fibers in the filename
  adjustFilename(filename, nFibersRow0);

  LOG(DEBUG) << "open file \"" << filename << "\".";
  // open file
  MPI_File fileHandle;
  MPIUtility::handleReturnValue(MPI_File_open(currentRankSubset_->mpiCommunicator(), filename.c_str(),
                                              //MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_UNIQUE_OPEN,
                                              MPI_MODE_WRONLY | MPI_MODE_CREATE,
                                              MPI_INFO_NULL, &fileHandle), "MPI_File_open");

  // write file header
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

  LOG(DEBUG) << "Fibers written to file \"" << filename << "\".";

  if (currentRankSubset_->ownRankNo() == 0)
  {
    fixInvalidFibersInFile(filename);

    if (nPointsWholeFiber != nNodesPerFiber_ || finalBottomZClip_ != bottomZClip_ || finalTopZClip_ != topZClip_)
    {
      resampleFibersInFile(nPointsWholeFiber, filename);
    }
  }
}

} // namespace
