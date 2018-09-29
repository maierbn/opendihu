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
    MPIUtility::handleReturnValue(MPI_File_write_ordered(fileHandle, writeBuffer.c_str(), writeBuffer.length(), MPI_BYTE, MPI_STATUS_IGNORE), "MPI_File_write_ordered");
  }
  else
  {
    MPIUtility::handleReturnValue(MPI_File_write_ordered(fileHandle, nullptr, 0, MPI_BYTE, MPI_STATUS_IGNORE), "MPI_File_write_ordered");
  }
}

};  // namespace
