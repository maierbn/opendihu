#include "gtest/gtest.h"

#include "easylogging++.h"

int argc;
char** argv;

namespace {
TEST(CompileTest, GTestCompiles)
{
  ASSERT_EQ(true, true);
}

}  // namespace

int main(int argc_, char **argv_) {
  ::testing::InitGoogleTest(&argc_, argv_);
  argc = argc_;
  argv = argv_;

  int returnValue = system("rm -rf SUCCESS1");
  if (returnValue != 0)
    LOG(INFO) << "failed to rm SUCCESS1 file.";

  int result = RUN_ALL_TESTS();

  // if all tests succeeded, create file "SUCCESS1"
  if (result == 0)
  {
    LOG(INFO) << "\x1b[32mAll tests passed on a single rank.\x1b[0m";
    std::ofstream outfile("SUCCESS1");
    outfile.close();
  }
  return result;
}
