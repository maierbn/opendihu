#include "postprocessing/parallel_fiber_estimation/parallel_fiber_estimation.h"

namespace Postprocessing
{

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
interpolateFineFibersFromFile()
{
  // open existing file to read
  bool waitIfFileGetsBig = specificSettings_.getOptionBool("waitIfFileGetsBig", true);

  struct stat statBuffer;
  int rc = stat(resultFilename_.c_str(), &statBuffer);
  if (rc != 0)
  {
    LOG(FATAL) << "Could not open file \"" << resultFilename_ << "\".";
  }
  int fileOldSize = statBuffer.st_size;

  std::ifstream fileOld(resultFilename_.c_str(), std::ios::in | std::ios::binary);
  assert (fileOld.is_open());

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
  memcpy(headerBuffer.data(), "opendihu fibers file            ", 32);

  int nFibersNewX = (nFibersOldX-1) * nFineGridFibers_ + nFibersOldX;
  int nFibersNew = MathUtility::sqr(nFibersNewX);

  std::stringstream newFilename;
  if (adjustFilename(resultFilename_, nFibersNewX))
  {
    newFilename << resultFilename_;
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

  if (nGibiBytes >= 1 && waitIfFileGetsBig)
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
      // index of bottom left fiber in old file
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
        /* if fiber is no key fiber, read 4 neighbouring fibers and interpolate
        // *     *
        //  \   /
        //    o
        //  /   \
        // x     *     // x = oldFiberIndex
        */
        // determine indices of 4 neighbouring fibers
        int oldFiberIndex[4] = {
          oldFiberIndexY * nFibersOldX + oldFiberIndexX,
          oldFiberIndexY * nFibersOldX + (oldFiberIndexX+1),
          (oldFiberIndexY+1) * nFibersOldX + oldFiberIndexX,
          (oldFiberIndexY+1) * nFibersOldX + (oldFiberIndexX+1)
        };

        // check if fiber is contained in file
        std::array<bool,4> edgeFiberIsValid({true,true,true,true});
        for (int i = 0; i < 4; i++)
        {
          // check if end of data for fiber is beyond end of file
          std::size_t oldFiberEndPosInFile = 32+headerLength + (oldFiberIndex[i]+1)*fiberDataSize;
          if (oldFiberEndPosInFile > fileOldSize)
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
        if (!edgeFiberIsValid[2] && !edgeFiberIsValid[3])
        {
          alpha1 = 0.0;
        }
        if (!edgeFiberIsValid[1] && !edgeFiberIsValid[3])
        {
          alpha0 = 0.0;
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
