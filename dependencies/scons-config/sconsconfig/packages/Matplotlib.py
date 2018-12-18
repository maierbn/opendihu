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
import matplotlib
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
# Matplotlib requires Python, Cython, NumpyC and bzip2 before Python
#
class Matplotlib(Package):
  
    def __init__(self, **kwargs):
        defaults = {
            'download_url': 'https://files.pythonhosted.org/packages/9e/59/f235ab21bbe7b7c6570c4abf17ffb893071f4fa3b9cf557b09b60359ad9a/matplotlib-2.2.3-cp36-cp36m-manylinux1_x86_64.whl'
        }
        defaults.update(kwargs)
        super(Matplotlib, self).__init__(**defaults)
        self.ext = '.cpp'
        self.sub_dirs = [
            ('include', 'lib'),
        ]
        self.headers = []
        self.libs = []
        self.check_text = check_text
        self.static = False
        
        if os.environ.get("PE_ENV") is not None:
          # Setup the build handler.
          self.set_build_handler([
              '$${DEPENDENCIES_DIR}/python/install/bin/pip3 install ${PREFIX}/../matplotlib-2.2.3-cp36-cp36m-manylinux1_x86_64.whl --prefix=${DEPENDENCIES_DIR}/python/install'
          ])
        else :
          # Setup the build handler.
          self.set_build_handler([
              '$${DEPENDENCIES_DIR}/python/install/bin/pip3 install matplotlib --prefix=${DEPENDENCIES_DIR}/python/install'
          ])
        
        # Matplotlib is installed in the directory tree of python, under lib/python3.6/site-packages. It is done using pip.

        self.number_output_lines = 13780
        
    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for Matplotlib ...    ')
        self.check_options(env)

        res = super(Matplotlib, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
