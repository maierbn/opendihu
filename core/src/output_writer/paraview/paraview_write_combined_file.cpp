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
    char b[1] = {' '};
    MPI_Status status;
    MPIUtility::handleReturnValue(MPI_File_write_ordered(fileHandle, b, 1, MPI_BYTE, &status), "MPI_File_write_ordered", &status);
  }
}

};  // namespace
