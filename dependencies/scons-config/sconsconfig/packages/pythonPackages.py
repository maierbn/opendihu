import sys, os
from distutils import sysconfig
from Package import Package

check_text = r'''
#include <stdlib.h>
#include <stdio.h>
#include <Python.h>
#include <iostream>

int main(int argc, char* argv[]) {

  Py_Initialize();
   return EXIT_SUCCESS;
}
'''

class pythonPackages(Package):
  
    def __init__(self, **kwargs):
        defaults = {
          'download_url': 'https://files.pythonhosted.org/packages/50/2f/618c5b6804e6dda90f024f0bc2d14ffc3db00221a818ee35da478427015d/svg.path-3.0-py2.py3-none-any.whl',
        }
        defaults.update(kwargs)
        super(pythonPackages, self).__init__(**defaults)
        self.ext = '.cpp'
        self.sub_dirs = [
            ('include', 'lib'),
        ]
        self.headers = []
        self.libs = []
        self.check_text = check_text
        self.static = False
        
        import socket
        
        if os.environ.get("PE_ENV") is not None:
          # Setup the build handler.
          self.set_build_handler([
              '$${DEPENDENCIES_DIR}/python/install/bin/pip3 install ${DEPENDENCIES_DIR}/pythonpackages/*.whl --prefix=${DEPENDENCIES_DIR}/python/install'
          ])
        elif socket.gethostname() == 'cmcs09':
           print("scons package pythonPackages.py: Nothing to be done.")
        else :
          # Setup the build handler.
          self.set_build_handler([
              '$${DEPENDENCIES_DIR}/python/install/bin/pip3 install numpy matplotlib scipy numpy-stl svg.path triangle --prefix=${DEPENDENCIES_DIR}/python/install'
          ])
        
        self.number_output_lines = 13780
        
    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for python packages...')
        self.check_options(env)

        res = super(pythonPackages, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
