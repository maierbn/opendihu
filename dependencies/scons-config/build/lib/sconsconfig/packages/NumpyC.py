import sys, os
from distutils import sysconfig
from Package import Package

check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <numpy/numpyconfig.h>

int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''


class NumpyC(Package):
  
    def __init__(self, **kwargs):
        defaults = {
            'download_url': 'https://github.com/numpy/numpy/archive/master.zip',
        }
        defaults.update(kwargs)
        super(NumpyC, self).__init__(**defaults)
        #self.ext = '.c'
        #self.sub_dirs = [
        #    ('include/mysql', 'lib'),
        #    ('include/mysql', 'lib64'),
        #]
        self.headers = ['numpy/numpyconfig.h']
        self.libs = [['libnpymath.a']]
        self.check_text = check_text
        self.static = False
        
        # Setup the build handler.
        self.set_build_handler([
            'export PYTHONPATH=$PYTHONPATH:${SOURCE_DIR}/../../../cython/install/lib/python2.7/site-packages/ && \
             python setup.py build_ext -i',
            'ln -s ${SOURCE_DIR}/numpy/core ${PREFIX}',
        ])

        self.number_output_lines = 2077
        
    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for Numpy C-API ... ')
        self.check_options(env)

        res = super(NumpyC, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
