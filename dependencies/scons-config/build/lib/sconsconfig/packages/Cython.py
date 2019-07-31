import sys, os
from distutils import sysconfig
from Package import Package

check_text = r'''
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[]) {
   return EXIT_SUCCESS;
}
'''

class Cython(Package):
  
    def __init__(self, **kwargs):
        defaults = {
            'download_url': 'https://pypi.python.org/packages/98/bb/cd2be435e28ee1206151793a528028e3dc9a787fe525049efb73637f52bb/Cython-0.27.2.tar.gz',
        }
        defaults.update(kwargs)
        super(Cython, self).__init__(**defaults)
        #self.ext = '.c'
        #self.sub_dirs = [
        #    ('include/mysql', 'lib'),
        #    ('include/mysql', 'lib64'),
        #]
        #self.headers = ['mysql.h']
        self.libs = []
        self.extra_libs = []
        self.check_text = check_text
        self.static = False
        
        # Setup the build handler.
        self.set_build_handler([
            'mkdir -p install',
            'export PYTHONPATH=$PYTHONPATH:${SOURCE_DIR}/install/lib/python2.7/site-packages/ && \
            python setup.py install --prefix=install',
        ])

        self.number_output_lines = 871
        
    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for Cython ... ')
        self.check_options(env)

        res = super(Cython, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
