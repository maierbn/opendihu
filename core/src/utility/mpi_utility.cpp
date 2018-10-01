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

    LOG(DEBUG) << "MPI status: count=" << status->_ucount << ", cancelled=" << status->_cancelled << ", MPI_SOURCE=" << status->MPI_SOURCE
      << ", MPI_TAG=" << status->MPI_TAG << ", MPI_ERROR=" << status->MPI_ERROR;
  }

  if (returnValue == MPI_SUCCESS)
    return;
  
  int stringLength = 200000;
  char errorString[stringLength];
  MPI_Error_string(returnValue, errorString, &stringLength);
  errorString[stringLength] = '\0';
  
  if (!descriptor.empty())
  {
    LOG(ERROR) << "Error in " << descriptor << ": " << std::string(errorString);
  }
  else 
  {
    LOG(ERROR) << "Error in MPI function: " << std::string(errorString);
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
    int pid = getpid();
    LOG(INFO) << "Rank " << rankNo << ", PID " << pid << " is waiting for gdbResume=" << gdbResume 
      << " to become 1 " << std::endl << std::endl
      << "sudo gdb -p " << pid << std::endl << std::endl
      << "select-frame 2" << std::endl
      << "set var gdbResume = 1" << std::endl
      << "info locals " << std::endl 
      << "continue";
    while (gdbResume == 0)
    {
      std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
    LOG(INFO) << "Rank " << rankNo << ", PID " << pid << " resumes because gdbResume=" << gdbResume;
  }
}

}  // namespace MPIUtility
