#pragma once

#include <Python.h>  // has to be the first included header
#include <mpi.h>
#include <iostream>

namespace MPIUtility
{

//! check if returnValue is an MPI_ERROR code, if yes print an error message, in descriptor the MPI functino name from which the returnValue was retrived, can be given.
void handleReturnValue(int returnValue, std::string descriptor, MPI_Status *status=MPI_STATUS_IGNORE);

//! Make all processes wait until one sets the local variable 'gdb_resume' to 1 from gdb
void gdbParallelDebuggingBarrier();

//! load the contents of the file in parallel, use case is when a short file should be read by a large number of ranks
std::string loadFile(std::string filename, MPI_Comm mpiCommunicator);

} // namespace MPIUtility
