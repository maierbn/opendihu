#include "gtest/gtest.h"
#include <mpi.h>

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

  if (ret == 0)
  {
    std::ofstream outfile("SUCCESS");
    outfile.close();
  }

  return nFails;
}
