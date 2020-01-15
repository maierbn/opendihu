import sys, os, multiprocessing, subprocess
from Package import Package

class bzip2(Package):

    def __init__(self, **kwargs):
        defaults = {
          'download_url': 'http://www.bzip.org/1.0.6/bzip2-1.0.6.tar.gz'
            
        }
        defaults.update(kwargs)
        super(bzip2, self).__init__(**defaults)
        self.ext = '.cpp'
        self.sub_dirs = [
            ('include', 'lib'),
        ]
        
        self.check_text = r'''
          #include <stdlib.h>
          #include <stdio.h>
          #include <bzlib.h>
          int main(int argc, char* argv[])
          {
            return EXIT_SUCCESS;
          }
        '''
    
        # Setup the build handler.
        self.libs = ["bz2"]
        self.headers = ["bzlib.h"]
        
        self.set_build_handler([
          'cd ${SOURCE_DIR} && make CFLAGS="-fPIC -Wall -Winline -O2 -g -D_FILE_OFFSET_BITS=64" && make install PREFIX=${PREFIX}'
        ])

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for bzip2 ... ')
        self.check_options(env)

        res = super(bzip2, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
