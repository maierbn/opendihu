#include "utility/mpi_utility.h"

#include <mpi.h>
#include "easylogging++.h"

#include <sys/types.h>  // getpid
#include <unistd.h>     // getpid
#include <thread>
#include <chrono>

namespace MPIUtility
{
  
void handleReturnValue(int returnValue, std::string descriptor, MPI_Status *status)
{
  if (status != MPI_STATUS_IGNORE && status != nullptr)
  {
    /*
     * typedef struct _MPI_Status {
  int count;
  int cancelled;
  int MPI_SOURCE;
  int MPI_TAG;
  int MPI_ERROR;
} MPI_Status, *PMPI_Status;
*/

    //LOG(DEBUG) << "MPI status: count=" << status->_ucount << ", cancelled=" << status->_cancelled << ", MPI_SOURCE=" << status->MPI_SOURCE
    //  << ", MPI_TAG=" << status->MPI_TAG << ", MPI_ERROR=" << status->MPI_ERROR;
  }

  if (returnValue == MPI_SUCCESS)
    return;
  
  int stringLength = 200000;
  std::vector<char> errorString(stringLength);
  MPI_Error_string(returnValue, errorString.data(), &stringLength);
  errorString[stringLength] = '\0';
  
  if (!descriptor.empty())
  {
    LOG(ERROR) << "Error in " << descriptor << ": " << std::string(errorString.begin(), errorString.begin()+stringLength);
  }
  else 
  {
    LOG(ERROR) << "Error in MPI function: " << std::string(errorString.begin(), errorString.begin()+stringLength);
  }
}
 
void gdbParallelDebuggingBarrier()
{
  volatile int gdbResume = 0;
  
  int nRanks, rankNo;
  handleReturnValue (MPI_Comm_size(MPI_COMM_WORLD, &nRanks));
  handleReturnValue (MPI_Comm_rank(MPI_COMM_WORLD, &rankNo));
  
  if (nRanks > 0) 
  {
#ifdef NDEBUG
    int pid = getpid();
    LOG(INFO) << "Rank " << rankNo << ", PID " << pid << " is waiting for gdbResume=" << gdbResume 
      << " to become 1 " << std::endl << std::endl
      << "gdb -p " << pid << std::endl << std::endl
      << "select-frame 2" << std::endl
      << "set var gdbResume = 1" << std::endl
      << "info locals " << std::endl 
      << "continue";
#else
    int pid = getpid();
    LOG(DEBUG) << "Rank " << rankNo << ", PID " << pid << " is waiting for gdbResume=" << gdbResume
      << " to become 1 " << std::endl << std::endl
      << "gdb -p " << pid << std::endl << std::endl
      << "select-frame 2" << std::endl
      << "set var gdbResume = 1" << std::endl
      << "info locals " << std::endl
      << "continue" << std::endl << std::endl;

#endif
    while (gdbResume == 0)
    {
      std::this_thread::sleep_for (std::chrono::milliseconds(5));
    }
    LOG(INFO) << "Rank " << rankNo << ", PID " << pid << " resumes because gdbResume=" << gdbResume;
  }
}

std::string loadFile(std::string filename, MPI_Comm mpiCommunicator)
{
  int nRanks, rankNo;
  handleReturnValue (MPI_Comm_size(mpiCommunicator, &nRanks));
  handleReturnValue (MPI_Comm_rank(mpiCommunicator, &rankNo));

  unsigned long long fileLength = 0;

  // on master rank, check if file exists and how many bytes it contains
  if (rankNo == 0)
  {
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    if (!in.is_open())
    {
      LOG(FATAL) << "Could not open file \""  << filename << "\".";
    }
    fileLength = in.tellg();
    in.close();
  }

  handleReturnValue(MPI_Bcast(&fileLength, 1, MPI_UNSIGNED_LONG_LONG, 0, mpiCommunicator), "MPI_Bcast");

  // collectively open the file for reading
  MPI_File fileHandle;
  handleReturnValue(MPI_File_open(mpiCommunicator, filename.c_str(),
                                  MPI_MODE_RDONLY,
                                  MPI_INFO_NULL, &fileHandle), "MPI_File_open");

  std::vector<char> fileContents(fileLength);
  int offset = 0;
  handleReturnValue(MPI_File_read_at_all(fileHandle, offset, fileContents.data(), fileLength, MPI_BYTE, MPI_STATUS_IGNORE), "MPI_Read_at_all");

  handleReturnValue(MPI_File_close(&fileHandle), "MPI_File_close");

  return std::string(fileContents.begin(), fileContents.end());
}

}  // namespace MPIUtility
