#include "postprocessing/scale_fibers_in_file.h"

#include "control/types.h"
#include "utility/math_utility.h"

namespace Postprocessing
{

void scaleFibersInFile(std::string inputFilename, std::string outputFilename, double scalingFactor)
{
  // create a new file with all the fibers from the old file but resampled such that they have nNodesPerFiber_ nodes
  LOG(DEBUG) << "scale all fibers in file, this is completely serial";

  // determine nFibers of the input file
  std::ifstream file(inputFilename.c_str(), std::ios::in | std::ios::binary);

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
  int nPointsPerFiber = bufferNPointsPerFiber.i;
  file.close();

  int nFibersX = int(std::round(std::sqrt(nFibers)));

  // open existing file to read
  std::ifstream fileOld(inputFilename.c_str(), std::ios::in | std::ios::binary);
  assert (fileOld.is_open());

  // open new file to write
  LOG(DEBUG) << "write to file " << outputFilename;
  std::ofstream fileNew(outputFilename.c_str(), std::ios::out | std::ios::binary);
  assert (fileNew.is_open());

  const int headerLength = sizeof(int32_t)*10;

  // copy header
  std::vector<char> headerBuffer(32+headerLength);
  fileOld.read(headerBuffer.data(), 32+headerLength);

  // set first 32 bytes
  strncpy(headerBuffer.data(), "opendihu fibers file            ", 32);
  fileNew.write(headerBuffer.data(), 32+headerLength);

  // write new date
  union
  {
    int32_t parameter;
    char c[sizeof(int32_t)];
  };
  parameter = time(NULL);

  fileNew.seekp(32+9*4);
  fileNew.write(c, 4);

  fileNew.seekp(32+headerLength);

  const long long fiberDataSize = nPointsPerFiber*3*sizeof(double);

  std::cout << "Scaling factor " << scalingFactor << ", " << nFibersX << "x" << nFibersX << " fibers with " << nPointsPerFiber << " points per fiber." << std::endl;

  Vec3 oldBoundingBoxMinimum;
  Vec3 oldBoundingBoxMaximum;
  Vec3 newBoundingBoxMinimum;
  Vec3 newBoundingBoxMaximum;

  // loop over fibers in old file
  for (int fiberIndex = 0; fiberIndex != nFibers; fiberIndex++)
  {
    // loop over nodes of the new fiber and write them to the new file
    for (int zIndex = 0; zIndex < nPointsPerFiber; zIndex++)
    {
      int fileIndex = 32+headerLength + fiberIndex*fiberDataSize + zIndex*3*sizeof(double);

      // read point from input file
      Vec3 oldPoint;
      fileOld.seekg(fileIndex);
      MathUtility::readPoint(fileOld, oldPoint);

      // compute new point
      Vec3 newPoint = oldPoint * scalingFactor;

      // write point to file
      MathUtility::writePoint(fileNew, newPoint);

      // update boundary box
      if (fiberIndex == 0 && zIndex == 0)
      {
        oldBoundingBoxMinimum = oldPoint;
        oldBoundingBoxMaximum = oldPoint;
        newBoundingBoxMinimum = newPoint;
        newBoundingBoxMaximum = newPoint;
      }
      for (int i = 0; i < 3; i++)
      {
        oldBoundingBoxMinimum[i] = std::min(oldBoundingBoxMinimum[i], oldPoint[i]);
        oldBoundingBoxMaximum[i] = std::max(oldBoundingBoxMaximum[i], oldPoint[i]);
        newBoundingBoxMinimum[i] = std::min(newBoundingBoxMinimum[i], newPoint[i]);
        newBoundingBoxMaximum[i] = std::max(newBoundingBoxMaximum[i], newPoint[i]);
      }
    }
  }

  fileOld.close();
  fileNew.close();

  std::cout << "Input file \"" << inputFilename << "\",\n  bounding box "
    << "[" << oldBoundingBoxMinimum[0] << "," << oldBoundingBoxMaximum[0] << "] x "
    << "[" << oldBoundingBoxMinimum[1] << "," << oldBoundingBoxMaximum[1] << "] x "
    << "[" << oldBoundingBoxMinimum[2] << "," << oldBoundingBoxMaximum[2] << "],\n"
    << "  size " << oldBoundingBoxMaximum[0] - oldBoundingBoxMinimum[0] << " x "
    << oldBoundingBoxMaximum[1] - oldBoundingBoxMinimum[1] << " x "
    << oldBoundingBoxMaximum[2] - oldBoundingBoxMinimum[2] << std::endl;

  std::cout << "Output file \"" << outputFilename << "\",\n  bounding box "
    << "[" << newBoundingBoxMinimum[0] << "," << newBoundingBoxMaximum[0] << "] x "
    << "[" << newBoundingBoxMinimum[1] << "," << newBoundingBoxMaximum[1] << "] x "
    << "[" << newBoundingBoxMinimum[2] << "," << newBoundingBoxMaximum[2] << "],\n"
    << "  size " << newBoundingBoxMaximum[0] - newBoundingBoxMinimum[0] << " x "
    << newBoundingBoxMaximum[1] - newBoundingBoxMinimum[1] << " x "
    << newBoundingBoxMaximum[2] - newBoundingBoxMinimum[2] << std::endl;

}

} // namespace
