#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
resampleFibersInFile(int nPointsPerFiber, std::string filename)
{
  // create a new file with all the fibers from the old file but resampled such that they have nNodesPerFiber_ nodes
  LOG(DEBUG) << "resample fibers in file, this is completely serial, nBoundaryPointsXNew_: " << nBoundaryPointsXNew_;
  LOG(INFO);  // newline

  // determine nFibersX of the input file
  std::ifstream file(filename.c_str(), std::ios::in | std::ios::binary);

  // skip first part of header
  file.seekg(32+4);
  union int32
  {
    char c[4];
    int32_t i;
  }
  bufferNFibers, bufferNPointsPerFiber;

  // get number of fibers
  file.read(bufferNFibers.c, 4);
  int nFibers = bufferNFibers.i;

  // get number of points per fiber
  file.read(bufferNPointsPerFiber.c, 4);
  if (nPointsPerFiber != bufferNPointsPerFiber.i)
  {
    LOG(FATAL) << "Invalid number of points per fiber: file \"" << filename << "\" contains " << bufferNPointsPerFiber.i
      << ", but " << nPointsPerFiber << " are required.";
  }
  file.close();

  int nFibersX = int(std::round(std::sqrt(nFibers)));
  std::string filenameExistingFile = filename;

  // adjust result filename to contain the correct number of fibers
  adjustFilename(filename, nFibersX);

  // if filename does not change, i.e. filename does not depend on number of fibers
  if (filename == filenameExistingFile)
  {
    filenameExistingFile = filename + std::string(".lowres");

    // rename existing file
    std::stringstream moveCommand;
    moveCommand << "mv " << filename << " " << filenameExistingFile;
    int ret = std::system(moveCommand.str().c_str());
    ret++;
    LOG(INFO) << nFibersX << "x" << nFibersX << " fibers with " << nPointsPerFiber << " points per fiber written to file \"" << filenameExistingFile << "\".";
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
  }

  // open existing file to read
  std::ifstream fileOld(filenameExistingFile.c_str(), std::ios::in | std::ios::binary);
  assert (fileOld.is_open());

  // open new file to write
  LOG(DEBUG) << "write to file " << filename;
  std::ofstream fileNew(filename.c_str(), std::ios::out | std::ios::binary);
  assert (fileNew.is_open());

  const int headerLength = sizeof(int32_t)*10;

  // copy header
  std::vector<char> headerBuffer(32+headerLength);
  fileOld.read(headerBuffer.data(), 32+headerLength);

  // set first 32 bytes
  strncpy(headerBuffer.data(), "opendihu fibers file            ", 32);
  fileNew.write(headerBuffer.data(), 32+headerLength);

  // write new number of points per fiber
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

  const long long fiberDataSize = nPointsPerFiber*3*sizeof(double);

  // loop over fibers in old file
  for (int fiberIndex = 0; fiberIndex != nFibers; fiberIndex++)
  {
    Vec3 previousPoint;
    Vec3 nextPoint;

    int oldZIndexPrevious0 = -1;    /// z index of the point that is currently loaded in previousPoint
    int oldZIndexNext0 = -1;        /// z index of the point that is currently loaded in nextPoint

    // loop over nodes of the new fiber and write them to the new file
    for (int zIndex = 0; zIndex < nNodesPerFiber_; zIndex++)
    {
      // compute the z value of the current point in the new fiber
      double currentZ = finalBottomZClip_ + zIndex * double(finalTopZClip_ - finalBottomZClip_) / (nNodesPerFiber_ - 1);

      //LOG(DEBUG) << "clip: " << finalBottomZClip_ << "," << finalTopZClip_ << " zIndex: " << zIndex << "/" << nNodesPerFiber_ << ", currentZ: " << currentZ;

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

      //if (fiberIndex < 10 || fiberIndex > nFibers-10)
      if (false)
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

  LOG(INFO) << nFibersX << "x" << nFibersX << " fibers with " << nNodesPerFiber_ << " points per fiber written to file \"" << filename << "\"." << std::endl;
  LOG(INFO) << "Use the following command to convert the file to various formats for inspection:\n  examine_bin_fibers.py \"" << filenameExistingFile << "\"" << std::endl;
}

} // namespace
