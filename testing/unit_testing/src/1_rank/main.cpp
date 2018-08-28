#include "gtest/gtest.h"

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
  return RUN_ALL_TESTS();
}
