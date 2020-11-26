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
import triangle
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

class pythonUtils(Package):
  
    def __init__(self, **kwargs):
        defaults = {
        }
        defaults.update(kwargs)
        super(pythonUtils, self).__init__(**defaults)
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
              '$${DEPENDENCIES_DIR}/python/install/bin/pip3 install ${PREFIX}/../*.whl --prefix=${DEPENDENCIES_DIR}/python/install'
          ])
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

        res = super(pythonUtils, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
