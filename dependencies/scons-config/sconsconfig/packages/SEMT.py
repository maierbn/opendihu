import sys, os, multiprocessing
from Package import Package

class SEMT(Package):

    def __init__(self, **kwargs):
        defaults = {
            #'download_url': 'https://github.com/st-gille/semt/archive/master.zip',
            'download_url': 'https://github.com/maierbn/semt/archive/master.zip',
        }
        defaults.update(kwargs)
        super(SEMT, self).__init__(**defaults)
        self.ext = '.cpp'
        self.sub_dirs = [
            ('include', 'lib'),
        ]
        #self.headers = ['mysql.h']
        #self.libs = ['mysqlclient']
        #self.extra_libs = ['lapack', 'blas']
        #self.build_flags = "-DSEMT_DISABLE_PRINT"
        self.check_text = r'''
#include <iostream>
using namespace std;

// Include namespace SEMT & global operators.
#include "semt/Semt.h"
// Include macros: INT, DINT, VAR, DVAR, PARAM, DPARAM
#include "semt/Shortcuts.h"
using namespace SEMT;

template<typename T>
std::ostream &operator<<(std::ostream &stream, const std::vector<T> &values)
{
  if (values.empty())
  {
    stream << "()";
    return stream;
  }

  stream << "(" << values[0];
  for (unsigned long i=1; i<values.size(); i++)
    stream << "," << values[i];
  stream << ")";
  return stream;
}

int main(int argc, char* argv[])
{
  DVAR(x2, 2);
  DINT(Two, 2);
  DPARAM(t1, 1)
  cout << (VAR(0) * VAR(1) - t1 + pow(x2, Two)) << endl;
  cout << deriv_t(pow(VAR(0) * x2, PARAM(0)), x2) << endl;
// output:
// (((x0 * x1) - t1) + (x2)^(2))
// (((x0 * x2))^(t0) * ((t0 * x0) / (x0 * x2)))
  return EXIT_SUCCESS;
}
'''
        # Setup the build handler.
        self.set_build_handler([
          "mkdir -p ${PREFIX}/include ",
          "sed -i.bak '150,166d' ${SOURCE_DIR}/semt/Semtfwd.h",    # remove bogus definition of operator<< for vectors
          "ln -s ${SOURCE_DIR}/semt ${PREFIX}/include/semt",
          "ln -s ${SOURCE_DIR}/loki ${PREFIX}/include/loki",
          #"cd ${PREFIX}/../build && cmake \
          #  -DBUILD_TYPE=Release \
          #  -DINSTALL_PREFIX=${PREFIX} \
          # -DLIBXML2_LIBRARIES=../../libxml2/install/lib/libxml2.a \
          #  -DLIBXML2_INCLUDE_DIR=../../libxml2/install/include \
          #  -DLIBCELLML_BUILD_SHARED=Off \
          #  -DLIBCELLML_COVERAGE=Off \
          #  -DLIBCELLML_MEMCHECK=Off \
          #  -DLIBCELLML_UNIT_TESTS=Off \
          #  ${SOURCE_DIR} && make && make install",
          #"ln -s ${PREFIX}/include/libcellml/module/libcellml ${PREFIX}/include/libcellml"
        ])
        self.number_output_lines = 84
          
        self.libs = []
        self.headers = ["semt/Semt.h", "semt/Shortcuts.h"]
        

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for SEMT ... ')
        self.check_options(env)

        res = super(SEMT, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
