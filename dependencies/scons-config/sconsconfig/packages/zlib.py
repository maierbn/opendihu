import sys, os, multiprocessing, subprocess
from .Package import Package

class zlib(Package):

    def __init__(self, **kwargs):
        defaults = {
          'download_url': 'https://zlib.net/zlib-1.2.11.tar.gz'
            
        }
        defaults.update(kwargs)
        super(zlib, self).__init__(**defaults)
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
        self.libs = ["z"]
        self.headers = ["zlib.h"]
        
        self.set_build_handler([
          'cd ${SOURCE_DIR} && ./configure --prefix=${PREFIX} && make install'
        ])

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for zlib  ... ')
        self.check_options(env)

        res = super(zlib, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
