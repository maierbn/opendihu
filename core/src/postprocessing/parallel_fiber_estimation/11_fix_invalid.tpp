#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

#include <cstdio>
namespace Postprocessing
{

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
fixInvalidFibersInFile(std::string filename)
{
  // open the file again and interpolate all missing fibers

  // copy existing file
  std::string filenameExistingFile = filename + std::string(".unfixed");
  std::stringstream moveCommand;
  moveCommand << "cp " << filename << " " << filenameExistingFile;
  int ret = std::system(moveCommand.str().c_str());
  ret++;
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));

  int nPointsPerFiber = 0;
  int nFibers = 0;
  int headerLength = 0;
  std::fstream file(filename.c_str(), std::ios::out | std::ios::in | std::ios::binary);
  if (!file.is_open())
  {
    LOG(FATAL) << "Could not open file \"" << filename << "\".";
  }
  else
  {
    // determine size of file
    struct stat statBuffer;
    stat(filename.c_str(), &statBuffer);
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

    LOG(DEBUG) << "headerLength: " << headerLength << ", nFibers: " << nFibers
      << ", nPointsPerFiber: " << nPointsPerFiber << ", nFibersX: " << nFibersX << ", fiberDataSize: " << fiberDataSize
      << ", fileSize: " << fileSize;

    if (int((fileSize-(32+headerLength)) / fiberDataSize) != nFibers)
    {
      LOG(ERROR) << "File \"" << filename << "\" states to have " << nFibers << " fibers in header, but actually has "
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
                    LOG(WARNING) << "Could not fix fiber (" << fiberIndexX << "," << fiberIndexY << ") / (" << nFibersX << "," << nFibersX << "),"
                     << " first point: " << firstPoint << ", neighbouring points: " << neighbouringFiberFirstPoint[0] << "," << neighbouringFiberFirstPoint[1];

                    int zIndex = int(nPointsPerFiber/2);
                    file.seekg(32+headerLength + validIndex0*fiberDataSize + zIndex*3*sizeof(double));
                    MathUtility::readPoint(file, neighbouringFiberFirstPoint[0]);

                    file.seekg(32+headerLength + validIndex1*fiberDataSize + zIndex*3*sizeof(double));
                    MathUtility::readPoint(file, neighbouringFiberFirstPoint[1]);
                  }

                  if (firstPoint[2] - neighbouringFiberFirstPoint[1][2] > 1e-12)
                  {
                    LOG(WARNING) << "Could not fix fiber (" << fiberIndexX << "," << fiberIndexY << ") / (" << nFibersX << "," << nFibersX << "),"
                     << " first point: " << firstPoint << ", neighbouring points: " << neighbouringFiberFirstPoint[0] << "," << neighbouringFiberFirstPoint[1];

                    int zIndex = int(nPointsPerFiber/2)+1;
                    file.seekg(32+headerLength + validIndex0*fiberDataSize + zIndex*3*sizeof(double));
                    MathUtility::readPoint(file, neighbouringFiberFirstPoint[0]);

                    file.seekg(32+headerLength + validIndex1*fiberDataSize + zIndex*3*sizeof(double));
                    MathUtility::readPoint(file, neighbouringFiberFirstPoint[1]);
                  }

                  if (firstPoint[2] - neighbouringFiberFirstPoint[1][2] > 1e-12)
                  {
                    LOG(ERROR) << "Could not fix fiber (" << fiberIndexX << "," << fiberIndexY << ") / (" << nFibersX << "," << nFibersX << "),"
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

    file.close();

    if (nFibersFixed == 0)
    {
      LOG(DEBUG) << "No invalid fibers were fixed, i.e. file did not change. Delete intermediate file " << filenameExistingFile << ", because it is has the same contents as " << filename;
      remove(filenameExistingFile.c_str());
    }
  }
}

} // namespace
