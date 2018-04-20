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
  
  std::string text = R"(
print("hello")

import sys
print("pythonpath:")
for p in sys.path:
  print(p)

import scipy
scipy.show_config()

import numpy as np
print(np.pi)

import pip
installed_packages = pip.get_installed_distributions()
installed_packages_list = sorted(["%s==%s" % (i.key, i.version) for i in installed_packages])
print("installed packages: ",installed_packages_list)

)";
  int ret = 0;
  try
  {
    std::cout << "run test program" << std::endl;
    ret = PyRun_SimpleString(text.c_str());
  }
  catch(...)
  {
    std::cout << "test program failed, ret=" << ret << std::endl;
  }
  std::cout << std::string(80, '-') << std::endl;
  
  // if there was an error in the python code
  if (ret != 0)
  {
    if (PyErr_Occurred())
    {
      // print error message and exit
      PyErr_Print();
      std::cout << "An error occurred in the python test script." << std::endl;
      exit(EXIT_FAILURE);
    }
    exit(EXIT_FAILURE);
  }
   return EXIT_SUCCESS;
}
'''


#
# SciPy requires Python and NumpyC and Cython
#
class SciPy(Package):
  
    def __init__(self, **kwargs):
        defaults = {
            'download_url': 'https://github.com/scipy/scipy/archive/master.zip'
        }
        defaults.update(kwargs)
        super(SciPy, self).__init__(**defaults)
        self.ext = '.cpp'
        self.sub_dirs = [
            ('include', 'lib'),
        #    ('include/mysql', 'lib64'),
        ]
        self.headers = []
        self.libs = []
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
             cd ${SOURCE_DIR} && \
             ${DEPENDENCIES_DIR}/python/install/bin/python3 setup.py build && \
             ${DEPENDENCIES_DIR}/python/install/bin/python3 setup.py install --prefix ${DEPENDENCIES_DIR}/python/install',
        ])
        
        # Scipy is installed in the directory tree of python, under lib/python3.6/site-packages. It does not create any .h or .a files.

        self.number_output_lines = 14142
        
    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for Scipy ... ')
        self.check_options(env)

        res = super(SciPy, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
