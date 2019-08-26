import sys, os, multiprocessing
from Package import Package

class Sphinx(Package):

    def __init__(self, **kwargs):
        defaults = {
            'download_url': 'https://github.com/sphinx-doc/sphinx/archive/master.tar.gz',
        }
        defaults.update(kwargs)
        super(Sphinx, self).__init__(**defaults)
        self.check_text = r'''
#include <stdlib.h>
#include <stdio.h>
int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''

        # get number of available processors
        p = multiprocessing.cpu_count()

        # Setup the build handler.
        self.set_build_handler([
        "pip install -t ${PREFIX} .",
          #"cd ${SOURCE_DIR} && ./configure --prefix=${PREFIX} --without-python && make && make install",
          #"mv ${PREFIX}/include/libxml2/* ${PREFIX}/include && rmdir ${PREFIX}/include/libxml2"
        ])
        self.number_output_lines = 44
        # creates /usr/local/bin/sphinx-build 

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for Sphinx ... ')
        self.check_options(env)

        res = super(Sphinx, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
