#include "gtest/gtest.h"
#include <mpi.h>
#include <petsc.h>

#include "easylogging++.h"
#include "../utility.h"

int argc;
char** argv;

int main(int argc_, char **argv_) {
  ::testing::InitGoogleTest(&argc_, argv_);
  argc = argc_;
  argv = argv_;

  int returnValue = system("rm -rf SUCCESS6");
  if (returnValue != 0)
    LOG(INFO) << "failed to rm SUCCESS6 file.";
    
  int ret = RUN_ALL_TESTS();
  LOG(INFO) << "return value of RUN_ALL_TESTS: " << ret << ", number of failed tests: " << nFails;

  // gather the return value of RUN_ALL_TESTS to process with rank 0
  int result = 1;
  PetscErrorCode ierr;
  ierr = MPI_Reduce(&ret, &result, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD); CHKERRQ(ierr);

  // if all tests succeeded an all ranks, create file "SUCCESS6"
  if (result == 0)
  {
    LOG(INFO) << "\x1b[32mAll tests passed on 6 ranks.\x1b[0m";
    std::ofstream outfile("SUCCESS6");
    outfile.close();
  }

  ierr = MPI_Finalize(); CHKERRQ(ierr);

  return nFails;
}
