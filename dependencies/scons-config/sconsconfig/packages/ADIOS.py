import sys, os, multiprocessing
from Package import Package
import subprocess

class ADIOS(Package):

  def __init__(self, **kwargs):
    defaults = {
        'download_url': 'https://github.com/ornladios/ADIOS2/archive/v2.3.0.zip',
    }
    defaults.update(kwargs)
    super(ADIOS, self).__init__(**defaults)
    self.ext = '.cpp'
    #self.sub_dirs = [
    #    ('include/mysql', 'lib'),
    #    ('include/mysql', 'lib64'),
    #]
    #self.headers = ['mysql.h']
    #self.libs = ['mysqlclient']
    #self.extra_libs = ['ADIOS', 'blas']
    #self.set_rpath = False    # do not use dynamic linkage
    self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <adios2.h>
int main(int argc, char* argv[]) {
return EXIT_SUCCESS;
}
'''

    # reference blas, cmake based, dynamic libraries
    self.set_build_handler([
      'mkdir -p ${PREFIX}',
      'cd ${SOURCE_DIR} && mkdir -p build && cd build && cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} -DCMAKE_BUILD_TYPE=RELEASE ..',
      'cd ${SOURCE_DIR}/build && make all install'
    ])
    
    self.number_output_lines = 4626
      
    self.libs = [["adios2_atl", "adios2_cmenet", "adios2_cmib", "adios2_cmmulticast", "adios2_cmselect", "adios2_cmsockets", "adios2_cmudp", "adios2_dill", "adios2_enet", "adios2_evpath", "adios2_ffs", "adios2_f", "adios2", "adios2_sst"]]
    self.headers = ["adios2.h"]

  def check(self, ctx):
    env = ctx.env
    ctx.Message('Checking for ADIOS ...         ')
    self.check_options(env)

    res = super(ADIOS, self).check(ctx)

    self.check_required(res[0], ctx)
    ctx.Result(res[0])
    return res[0]
