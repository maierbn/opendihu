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


#
# NumpyC requires Python and Cython
#
class NumpyC(Package):
  
    def __init__(self, **kwargs):
        defaults = {
            #'download_url': 'https://github.com/numpy/numpy/archive/master.zip',
            'download_url': 'https://files.pythonhosted.org/packages/94/b8/09db804ddf3bb7b50767544ec8e559695b152cedd64830040a0f31d6aeda/numpy-1.14.4.zip',
        }
        defaults.update(kwargs)
        super(NumpyC, self).__init__(**defaults)
        #self.ext = '.c'
        self.sub_dirs = [
            ('include', 'lib'),
        #    ('include/mysql', 'lib64'),
        ]
        self.headers = ['numpy/numpyconfig.h']
        self.libs = [['libnpymath.a']]
        self.check_text = check_text
        self.static = False
        
        # Setup the build handler.
        self.set_build_handler([
            'cd ${SOURCE_DIR} && echo "[openblas]" > site.cfg',
            'cd ${SOURCE_DIR} && echo "libraries = openblas" >> site.cfg',
            'cd ${SOURCE_DIR} && echo "library_dirs = ${DEPENDENCIES_DIR}/lapack/install/lib" >> site.cfg',
            'cd ${SOURCE_DIR} && echo "include_dirs = ${DEPENDENCIES_DIR}/lapack/install/include" >> site.cfg',
            'cd ${SOURCE_DIR} && echo "runtime_library_dirs = ${DEPENDENCIES_DIR}/lapack/install/lib" >> site.cfg',
            '$export PYTHONPATH=$PYTHONPATH:${DEPENDENCIES_DIR}/python/install/lib/$(basename $(find ${DEPENDENCIES_DIR}/python/install/lib/ -maxdepth 1 -type d -name "python*"))/site-packages/ && \
             export PATH=${DEPENDENCIES_DIR}/cython/install/bin:$PATH && \
             cd ${SOURCE_DIR} && \
             ${DEPENDENCIES_DIR}/python/install/bin/python3 setup.py build && \
             ${DEPENDENCIES_DIR}/python/install/bin/python3 setup.py install --prefix ${DEPENDENCIES_DIR}/python/install',
            '$cd ${DEPENDENCIES_DIR}/python/install/lib/ \
             && cd $(find . -maxdepth 1 -type d -name "python*") \
             && cd site-packages/ \
             && cd $(find . -maxdepth 1 -type d -name "numpy*.egg") \
             && cd numpy/core \
             && rm -f ${PREFIX} && ln -s $(pwd) ${PREFIX}'
            #'ln -s ${SOURCE_DIR}/numpy/core ${PREFIX}',
        ])
        
        # Numpy is installed in the directory tree of python, under lib/python3.6/site-packages
        # The headers and libraries are also positioned there and then linked from the numpyc/install directory.

        self.number_output_lines = 2960
        
    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for Numpy C-API ... ')
        self.check_options(env)

        res = super(NumpyC, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
