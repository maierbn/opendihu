import sys, os, multiprocessing
from Package import Package

##
##  Handles LAPACK library and also BLAS. There are fortran and C versions for each.
## 
## - ``libblas-dev``: reference BLAS (not very optimized)
## - ``libatlas-base-dev``: generic tuned ATLAS, it is recommended to tune it to
##   the available hardware, see /usr/share/doc/libatlas3-base/README.Debian for
##   instructions
## - ``libopenblas-base``: fast and runtime detected so no tuning required but a
##   very recent version is needed (>=0.2.15 is recommended).  Older versions of
##   OpenBLAS suffered from correctness issues on some CPUs.
##
##   We use OpenBLAS
## 
##
class LAPACK(Package):

    def __init__(self, **kwargs):
        defaults = {
            #'download_url': 'https://github.com/Reference-LAPACK/lapack-release/archive/lapack-3.7.1.zip',
            'download_url': 'http://github.com/xianyi/OpenBLAS/archive/v0.2.20.tar.gz',
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
        #self.set_rpath = False    # do not use dynamic linkage
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

        use_reference_blas = False
        # Setup the build handler.
        
        
        if os.environ.get("LIBSCI_BASE_DIR") is not None:
          self.libs = ["sci_cray_mpi_mp"]
          print("Cray environment detected, using \"sci_cray_mpi_mp\" for LAPACK")

        elif False:
          # reference blas, make based, only static libraries
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
          
        elif use_reference_blas:
            
          # reference blas, cmake based, dynamic libraries
          self.set_build_handler([
            'mkdir -p ${PREFIX}',
            'cd ${PREFIX} && cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_SHARED_LIBS=ON -DCBLAS=ON -DLAPACKE=ON -DBUILD_TESTING=OFF ${SOURCE_DIR}',
            'cd ${PREFIX} && make',
            'mkdir -p ${PREFIX}/include && cp ${SOURCE_DIR}/*/*/*.h ${PREFIX}/include',
          ])
          
          self.number_output_lines = 4626
            
          self.libs = ["lapack", "blas"]
          self.headers = ["lapacke.h"]
        
        else:  
          # use OpenBLAS
          self.set_build_handler([
            'mkdir -p ${PREFIX}',
            'cd ${SOURCE_DIR} && make && make install PREFIX=${PREFIX}',
          ])
          self.number_output_lines = 18788
        
          self.libs = ["openblas"]
          self.headers = ["cblas.h"]
          self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <cblas.h>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for LAPACK ... ')
        self.check_options(env)

        res = super(LAPACK, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
