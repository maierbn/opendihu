import sys, os, multiprocessing
from Package import Package

class googletest(Package):

    def __init__(self, **kwargs):
        defaults = {
            'download_url': 'https://github.com/google/googletest/archive/master.zip',
        }
        defaults.update(kwargs)
        super(googletest, self).__init__(**defaults)
        self.ext = '.cpp'
        #self.sub_dirs = [
        #    ('include/mysql', 'lib'),
        #    ('include/mysql', 'lib64'),
        #]
        #self.headers = ['mysql.h']
        #self.libs = ['mysqlclient']
        #self.extra_libs = ['lapack', 'blas']
        self.check_text = r'''
#include "gtest/gtest.h"

namespace {
// Tests that the Foo::Bar() method does Abc.
TEST(CompileTest, GTestCompiles) {
  ASSERT_EQ(true, true);
}

}  // namespace

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
'''

        # get number of available processors
        p = multiprocessing.cpu_count()

        # Setup the build handler.
         
        self.set_build_handler([
          'cd ${SOURCE_DIR}/googletest/make && make && make gtest.a',
          'cd ${SOURCE_DIR}/googletest/make && ./sample1_unittest',
          'mkdir -p ${PREFIX}/install',
          'mkdir -p ${PREFIX}/lib',
          'ln -s ${SOURCE_DIR}/googletest/include ${PREFIX}/include',
          'ln -s ${SOURCE_DIR}/googletest/make/gtest.a ${PREFIX}/lib/libgtest.a',
        ])
          
        self.libs = ["gtest"]
        self.headers = ["gtest/gtest.h"]

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for googletest ... ')
        self.check_options(env)

        res = super(googletest, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
