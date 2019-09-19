import sys, os, multiprocessing, subprocess
from Package import Package

class MegaMol(Package):

    def __init__(self, **kwargs):
        defaults = {
          'download_url': 'https://github.com/UniStuttgart-VISUS/megamol/archive/master.zip'
        }
        defaults.update(kwargs)
        super(MegaMol, self).__init__(**defaults)
        self.ext = '.cpp'
        self.sub_dirs = [
            ('include', 'lib'),
        ]
        #self.headers = ['mysql.h']
        #self.libs = ['mysqlclient']
        #self.extra_libs = ['lapack', 'blas']
        self.check_text = r'''
          #include <stdlib.h>
          #include <stdio.h>
          //#include <MegaMol.h>
          int main(int argc, char* argv[])
          {
          
            return EXIT_SUCCESS;
          }
'''
    
        # Setup the build handler.
        #self.libs = ["MegaMol3.6m"]
        #self.headers = ["mmcore/Module.h"]
        
        # check
        self.set_build_handler([
          'mkdir -p ${PREFIX}',
          'cd ${SOURCE_DIR} && mkdir -p build && cd build && cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} .. && make && make install',
        ])
          
        self.number_output_lines = 7082

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for MegaMol ...       ')
        self.check_options(env)

        res = super(MegaMol, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
