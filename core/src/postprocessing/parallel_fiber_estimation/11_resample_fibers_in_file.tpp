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

template<typename BasisFunctionType>
void ParallelFiberEstimation<BasisFunctionType>::
stretchMeshAtCornersInFile(std::string filename)
{
  std::stringstream command;
  command << "cp " << filename << " " << filename << ".before_corner_stretch";
  int returnValue = system(command.str().c_str());
  if (returnValue != 0)
    LOG(DEBUG) << "Could not execute command \"" << command.str() << "\".";

  // open existing file to read
  std::fstream file(filename.c_str(), std::ios::in | std::ios::out | std::ios::binary);
  assert (file.is_open());

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
  int nFibersX = int(std::round(std::sqrt(nFibers)));

  // get number of points per fiber
  file.read(bufferNPointsPerFiber.c, 4);
  int nPointsPerFiber = bufferNPointsPerFiber.i;

  const int headerLength = sizeof(int32_t)*10;
  const long long fiberDataSize = nPointsPerFiber*3*sizeof(double);

  std::vector<Vec3> pointsOnPlane(nFibersX*nFibersX);

  // loop over nodes in z direction
  for (int zIndex = 0; zIndex < nPointsPerFiber; zIndex++)
  {
    Vec2 center{0,0};

    // read all points
    // loop over fibers in file
    for (int fiberIndexJ = 0; fiberIndexJ != nFibersX; fiberIndexJ++)
    {
      for (int fiberIndexI = 0; fiberIndexI != nFibersX; fiberIndexI++)
      {
        int fiberIndex = fiberIndexJ*nFibersX + fiberIndexI;
        Vec3 point;

        file.seekg(32+headerLength + fiberIndex*fiberDataSize + zIndex*3*sizeof(double));
        MathUtility::readPoint(file, point);

        pointsOnPlane[fiberIndex] = point;

        center[0] += point[0];
        center[1] += point[1];
      }
    }
    center[0] /= nFibersX*nFibersX;
    center[1] /= nFibersX*nFibersX;

    // compute transformation
    // loop over fibers in file
    for (int fiberIndexJ = 0; fiberIndexJ != nFibersX; fiberIndexJ++)
    {
      for (int fiberIndexI = 0; fiberIndexI != nFibersX; fiberIndexI++)
      {
        // computations in index space
        double indexCenter = (nFibersX-1)/2;
        double indexI = fiberIndexI - indexCenter;
        double indexJ = fiberIndexJ - indexCenter;
        double indexPhi = atan2(indexJ, indexI);
        double a1 = std::fmod(indexPhi+2*M_PI, M_PI_2) / M_PI_2;

        // the following computations calculate two values: a in [0,1] and factor around 1
        // a specifies how much the factor is used and how much the normal radius is used, a = 0 means normal radius, i.e., no mesh change, a=1 means full changed mesh

        double h = 0.6;   // parameter, how far the area of mesh changes reaches in circumferential direction, 0<h<=1, smaller is more concentrated around corners
        double a = (1-sin(2*M_PI*((a1-0.5)/h+0.5)+M_PI_2))/2;    // 0 -> 0  0.5 -> 1  1 -> 0

        if (a1 < 0.5-h/2 || a1 > 0.5+h/2)
          a = 0;

        double indexR = sqrt(indexI*indexI + indexJ*indexJ);

        double maxR = fabs((nFibersX-1 - indexCenter)/cos(indexPhi));  // cos(phi) = a/h = x/maxR => maxR = x/cos(phi)
        if ((indexPhi >= M_PI_4 && indexPhi < 3*M_PI_4) || (indexPhi <= -M_PI_4 && indexPhi > -3*M_PI_4))
          maxR = fabs((nFibersX-1 - indexCenter)/sin(indexPhi));  // sin(phi) = g/h = y/maxR => maxR = y/sin(phi)

        double relativeR = fabs(indexR / maxR);
        double m = 0.9;   // 0.9        // parameter 0 << m <= 1, how much the edge elements will be stretches, higher m means more
        double xx = (relativeR*m);
        double scalingFunction1 = asin(xx) / asin(m);

        double scalingFunction2 = relativeR;

        double s = 0.5;       // [1-s, 1]    // parameter 0 < s << 1, how far the area of mesh changes reaches in radial direction, smaller is more concentrated around the corners
        double alpha = 0.7;   // 0 << alpha < 1, parameter how much the mesh is pulled towards the center, lower value = more, 1 = not at all

        double s2 = std::pow(s, 2);
        if (relativeR > 1-s)
        {
          // the following function is derived in doc/sympy/extend_mesh.py
          double x = relativeR;
          scalingFunction2 = 4*alpha - 8*alpha/s + 4*alpha/s2 + x*(-4*alpha + 16*alpha/s - 12*alpha/s2 + 5 - 16/s + 12/s2) - 4 + 8/s + std::pow(x, 3)*(4 - 4*alpha)/s2 + std::pow(x, 2)*(-8*alpha*s + 12*alpha + 8*s - 12)/s2 - 4/s2;
        }

        //double scalingFunction = 0.5*(scalingFunction1 + scalingFunction2);
        double scalingFunction = scalingFunction2;

        if (relativeR < 1e-10)
        {
          scalingFunction = 0;
          relativeR = 1;
        }

        double factor = scalingFunction / relativeR;  // map [0,1] -> [0,1]

        // do not change the border points
        if (fiberIndexI == 0 || fiberIndexI == nFibersX-1 || fiberIndexJ == 0 || fiberIndexJ == nFibersX-1)
        {
          a = 0;
          factor = 1;
        }

        // computations in world space
        int fiberIndex = fiberIndexJ*nFibersX + fiberIndexI;
        Vec3 point = pointsOnPlane[fiberIndex];

        // compute cylindrical coordinates
        double x = point[0] - center[0];
        double y = point[1] - center[1];
        double r = sqrt(x*x + y*y);
        double phi = atan2(y,x);

        // compute the transformed radius
        double newR = (1-a)*r + a*factor*r;

        pointsOnPlane[fiberIndex][0] = center[0] + newR*cos(phi);
        pointsOnPlane[fiberIndex][1] = center[1] + newR*sin(phi);
      }
    }

    // write points back to file
    // loop over fibers in file
    for (int fiberIndexJ = 0; fiberIndexJ != nFibersX; fiberIndexJ++)
    {
      for (int fiberIndexI = 0; fiberIndexI != nFibersX; fiberIndexI++)
      {
        int fiberIndex = fiberIndexJ*nFibersX + fiberIndexI;

        // write point to file
        file.seekg(32+headerLength + fiberIndex*fiberDataSize + zIndex*3*sizeof(double));
        MathUtility::writePoint(file, pointsOnPlane[fiberIndex]);
      }
    }
  }

  file.close();
}

} // namespace
