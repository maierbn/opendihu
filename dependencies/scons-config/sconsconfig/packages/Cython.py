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
            'download_url': 'https://files.pythonhosted.org/packages/b3/ae/971d3b936a7ad10e65cb7672356cff156000c5132cf406cb0f4d7a980fd3/Cython-0.28.3.tar.gz'
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

        self.number_output_lines = 883
        
    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for Cython ... ')
        self.check_options(env)

        res = super(Cython, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
