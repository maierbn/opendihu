import sys, os, multiprocessing
from Package import Package

##
##  Handles LAPACK library and also BLAS. There are fortran and C versions for each.
##
class LAPACK(Package):

    def __init__(self, **kwargs):
        defaults = {
            'download_url': 'https://github.com/Reference-LAPACK/lapack-release/archive/lapack-3.7.1.zip',
        }
        defaults.update(kwargs)
        super(LAPACK, self).__init__(**defaults)
        #self.ext = '.c'
        #self.sub_dirs = [
        #    ('include/mysql', 'lib'),
        #    ('include/mysql', 'lib64'),
        #]
        #self.headers = ['mysql.h']
        #self.libs = ['mysqlclient']
        #self.extra_libs = ['lapack', 'blas']
        self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <lapacke.h>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''

        # get number of available processors
        p = multiprocessing.cpu_count()

        # Setup the build handler.
        
        # make based, only static libraries
        if False:
          self.set_build_handler([
              'cp make.inc.example make.inc',
              'make lapack_install -j '+str(p),
              'make lapacklib -j '+str(p),
              'make tmglib -j '+str(p),
              'make blaslib -j '+str(p),
              'make cblaslib -j '+str(p),
              'make lapackelib -j '+str(p),
              'mkdir ${PREFIX}',
              'ln -s ${SOURCE_DIR}/LAPACKE/include/ ${PREFIX}/include',
              'mkdir ${PREFIX}/lib',
              'ln -s ${SOURCE_DIR}/libcblas.a ${PREFIX}/lib/',
              'ln -s ${SOURCE_DIR}/liblapack.a ${PREFIX}/lib/',
              'ln -s ${SOURCE_DIR}/liblapacke.a ${PREFIX}/lib/',
              'ln -s ${SOURCE_DIR}/librefblas.a ${PREFIX}/lib/',
              'ln -s ${PREFIX}/lib/librefblas.a ${PREFIX}/lib/libblas.a',
              'ln -s ${SOURCE_DIR}/libtmglib.a ${PREFIX}/lib/',
          ])
          self.number_output_lines = 4768
          
        # cmake based, dynamic libraries
        self.set_build_handler([
          'mkdir ${PREFIX}',
          'cd ${PREFIX} && cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_SHARED_LIBS=ON -DLAPACKE=ON ${SOURCE_DIR}',
          '!cd ${PREFIX} && make && make install',
          'cd ${PREFIX} && make && make install',
          'mkdir ${PREFIX}/include && cp ${SOURCE_DIR}/*/*/*.h ${PREFIX}/include',
        ])
        self.number_output_lines = 3399
        
        self.libs = ["lapack", "blas"]
        self.headers = ["lapacke.h"]

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for LAPACK ... ')
        self.check_options(env)

        res = super(LAPACK, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
