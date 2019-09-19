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
import stl
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

class NumpySTL(Package):
  
    def __init__(self, **kwargs):
        defaults = {
            'download_url': 'https://files.pythonhosted.org/packages/62/de/8362b69f9e0a13dfffd62a12e511bdeb2aa91c2d53b80017446eb561d3e2/numpy_stl-2.3.2-py3-none-any.whl'
        }
        defaults.update(kwargs)
        super(NumpySTL, self).__init__(**defaults)
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
              '$${DEPENDENCIES_DIR}/python/install/bin/pip3 install ${PREFIX}/../numpy_stl-2.3.2-py3-none-any.whl --prefix=${DEPENDENCIES_DIR}/python/install'
          ])
        else :
          # Setup the build handler.
          self.set_build_handler([
              '$${DEPENDENCIES_DIR}/python/install/bin/pip3 install numpy-stl --prefix=${DEPENDENCIES_DIR}/python/install'
          ])
        
        # NumpySTL is installed in the directory tree of python, under lib/python3.6/site-packages. It is done using pip.

        self.number_output_lines = 13780
        
    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for Numpy-stl ...     ')
        self.check_options(env)

        res = super(NumpySTL, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
