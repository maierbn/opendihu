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
  int ret = RUN_ALL_TESTS();
  LOG(INFO) << "return value of RUN_ALL_TESTS: " << ret << ", number of failed tests: " << nFails;

  // gather the return value of RUN_ALL_TESTS to process with rank 0
  int result = 1;
  PetscErrorCode ierr;
  ierr = MPI_Reduce(&ret, &result, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD); CHKERRQ(ierr);

  // if all tests succeeded an all ranks, create file "SUCCESS"
  if (result == 0)
  {
    LOG(INFO) << "All tests passed on all ranks.";
    std::ofstream outfile("SUCCESS");
    outfile.close();
  }

  ierr = MPI_Barrier(MPI_COMM_WORLD); CHKERRQ(ierr);

  return nFails;
}
