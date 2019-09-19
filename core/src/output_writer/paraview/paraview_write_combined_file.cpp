#include "output_writer/paraview/paraview.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

#include "easylogging++.h"
#include "base64.h"

#include "output_writer/paraview/loop_output.h"
#include "output_writer/paraview/loop_collect_mesh_properties.h"
#include "output_writer/paraview/poly_data_properties_for_mesh.h"

namespace OutputWriter
{

void Paraview::writeAsciiDataShared(MPI_File fileHandle, int ownRankNo, std::string writeBuffer)
{
  // collective blocking write, only rank 0 writes, but afterwards all have the same shared file pointer position
  if (ownRankNo == 0)
  {
    MPI_Status status;
    MPIUtility::handleReturnValue(MPI_File_write_ordered(fileHandle, writeBuffer.c_str(), writeBuffer.length(), MPI_BYTE, &status), "MPI_File_write_ordered", &status);
  }
  else
  {
    char b[1];
    MPI_Status status;
    MPIUtility::handleReturnValue(MPI_File_write_ordered(fileHandle, b, 0, MPI_BYTE, &status), "MPI_File_write_ordered", &status);
  }
}

void Paraview::writeCombinedTypesVector(MPI_File fileHandle, int ownRankNo, int nValues, bool output3DMeshes, int identifier)
{
  std::string writeBuffer;

  if (binaryOutput_)
  {
    if (output3DMeshes)
    {
      std::vector<int> values(nValues, 12);
      writeBuffer = Paraview::encodeBase64UInt8(values.begin(), values.end());
    }
    else
    {
      std::vector<int> values(nValues, 9);
      writeBuffer = Paraview::encodeBase64UInt8(values.begin(), values.end());
    }
  }
  else
  {
    for (int i = 0; i < nValues; i++)
    {
      if (output3DMeshes)
      {
        writeBuffer += std::string("12 ");
      }
      else
      {
        writeBuffer += std::string("9 ");
      }
    }
  }

  // collective blocking write, only rank 0 writes, but afterwards all have the same shared file pointer position
  if (ownRankNo == 0)
  {
    MPI_Status status;
    MPIUtility::handleReturnValue(MPI_File_write_ordered(fileHandle, writeBuffer.c_str(), writeBuffer.length(), MPI_BYTE, &status), "MPI_File_write_ordered", &status);
  }
  else
  {
    char b[1];
    MPI_Status status;
    MPIUtility::handleReturnValue(MPI_File_write_ordered(fileHandle, b, 0, MPI_BYTE, &status), "MPI_File_write_ordered", &status);
  }
}

//! constructor, initialize nPoints and nCells to 0
Paraview::VTKPiece::VTKPiece()
{
  properties.nPointsLocal = 0;
  properties.nCellsLocal = 0;
  properties.nPointsGlobal = 0;
  properties.nCellsGlobal = 0;
  properties.dimensionality = 0;
}

//! assign the correct values to firstScalarName and firstVectorName, only if properties has been set
void Paraview::VTKPiece::setVTKValues()
{
  // set values for firstScalarName and firstVectorName from the values in pointDataArrays
  for (auto pointDataArray : properties.pointDataArrays)
  {
    if (firstScalarName == "" && pointDataArray.second == 1)
      firstScalarName = pointDataArray.first;
    if (firstVectorName == "" && pointDataArray.second != 1)
      firstVectorName = pointDataArray.first;
  }
}

} // namespace
