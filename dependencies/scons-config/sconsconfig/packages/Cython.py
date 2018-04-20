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
        self.sub_dirs = [
            ('include', 'lib'),
        ]
        #self.headers = ['mysql.h']
        self.libs = []
        self.extra_libs = []
        self.check_text = check_text
        self.static = False
        
        # Setup the build handler.
        self.set_build_handler([
            'mkdir -p ${PREFIX}/include',
            '$export PYTHONPATH=$PYTHONPATH:${DEPENDENCIES_DIR}/python/install/lib/$(basename $(find ${DEPENDENCIES_DIR}/python/install/lib/ -maxdepth 1 -type d -name "python*"))/site-packages/ && \
            ${DEPENDENCIES_DIR}/python/install/bin/python3 setup.py install --prefix ${DEPENDENCIES_DIR}/python/install',
            'ln -s ${SOURCE_DIR}/bin ${PREFIX}/bin',
            '$export PATH=${PREFIX}/bin:$PATH'
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
